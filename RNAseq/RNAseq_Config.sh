#!/bin/bash
##//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## By: Jessica Qiu
## Supervisor: Dr. Rooksana Noorai, Dr. Vijay Shankar, and Dr. Christopher Saski
## Most Current Update Date: September 14, 2017
## Description: The purpose of this script is to allow the user to type in their inputs from the command line. The  script will create
##              multiple jobs to run all of the RNA sequence samples that utilize bioinformatic software such as FastQC, Trimmomatic, GSnap,
##              and Subread Feature counts.  More information on each tool can be found in RNAseq_UI-SampleFastQCToFeatureCounts.pbs.
##
## Pre-Setup:   Make sure the subfolder containing your .fastq or .fastq.gz files is called "raw_data".
## How To Run:
##              1. Place the RNAseq_UI-SampleFastQCToFeatureCounts.pbs. and the RNAseq_UI.sh in your script folder.
##              2. Adjust the configuration file, "RNAseqsAuto.config", to your information/inputs
##              3. Run RNAseq_UI.sh using the command "./RNAseq_UI" or "bash RNAseq_UI.sh".
##
## Close: Thank you to the great staff at the Clemson University Genomics Institute and the Clemson University Professional Internship/
##        Co-op Program for providing me with the opportunity to learn more about bioinformatics! And special thanks to Chelky Lin!
##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

echo "Total Number of Configuration Files Inputted: $#"

MODULES_NEEDED=$(cat <<MOD_NEED
module add java/1.8.0
module add fastqc/0.11.5
module add samtools/1.3.1
module add gcc/6.3.0
module add zlib/1.2.8
module add mpich/3.1.4
module add bowtie2/2.2.8
module add Subread/1.5.3
module add gmap_gsnap
module add perl/5.26.0
MOD_NEED
)

FASTQC_1=$(cat <<FASTQC
echo "Running FastQC..."
fastqc \${filename}*.fastq.gz

ls -lrth

time cp *.html \${workingdir}/1_FastQC
time cp *.zip \${workingdir}/1_FastQC
FASTQC
)

TRIM_2=$(cat <<TRIM
trimjar="/zfs/gcl/software/trimmomatic/0.36/trimmomatic-0.36.jar"
illumina_clip="/zfs/gcl/software/trimmomatic/0.36/adapters/illuminaClipping.fa"

time cp \${trimjar} \${illumina_clip} \${localdir1}
echo "Trimming for paired-end data"
time java -jar trimmomatic-0.36.jar PE -threads 16 -phred33\
 \${filename}_R1_001.fastq.gz\
 \${filename}_R2_001.fastq.gz\
 \${filename}_R1.paired.fastq.gz\
 \${filename}_R1.unpaired.fastq.gz\
 \${filename}_R2.paired.fastq.gz\
 \${filename}_R2.unpaired.fastq.gz\
 ILLUMINACLIP:illuminaClipping.fa:2:30:10\
 LEADING:3 TRAILING:3 HEADCROP:3 SLIDINGWINDOW:4:15 MINLEN:36

ls -lrth

## Copy trimmed fastq.gz files
echo "Copying trimmed fastq.gz files..."

time cp \${filename}_R1.paired.fastq.gz \${workingdir}/2_Trim
time cp \${filename}_R2.paired.fastq.gz \${workingdir}/2_Trim

## STEP 2.5: FASTQC ON TRIMMED DATA
echo "Running FastQC for the trimmed data..."
fastqc \${filename_forward}
fastqc \${filename_reverse}

time cp *.html \${workingdir}/2.5_PostTrimFastQC
time cp *.zip \${workingdir}/2.5_PostTrimFastQC

ls -lrth
TRIM
)

GSNAP_3=$(cat <<GSNAP
echo "START GSNAP"
cp \${workingdir}/2_Trim/* \${localdir1}

PATH=/zfs/gcl/software/gmap-gsnap/2017-11-15/bin:$PATH

## Align to the reference genome
echo "Aligning to the reference genome..."
time gsnap --gunzip -d \${refname}\
 \${filename_forward} \${filename_reverse}\
 -t 15\
 --orientation=RF\
 -N 1\
 -A sam\
 > \${filename}.sam

echo "Removing fastq file..."
rm \${filename_forward}
rm \${filename_reverse}

## Converting sam to bam...
echo "Converting sam to bam file..."
time samtools view -bS\
 \${filename}.sam\
 > \${filename}.bam
echo "Conversion Completed."
ls -lrth

rm \${filename}.sam

## Sorting the bamfile
echo "Sorting the bamfile..."
time samtools sort\
 \${filename}.bam\
 -o \${filename}.sorted.bam
echo "Sorting Completed."

## Indexing the bamfile
echo "Indexing the bamfile..."
time samtools\
 index\
 \${filename}.sorted.bam
echo "Indexing Completed"
ls -lrth

rm \${filename}.bam

## Copying output files back to /zfs/gcl
echo "Copying output files back to /zfs/gcl..."
time cp \${filename}* \${resultsdir}
GSNAP
)

FEATURECOUNTS_4=$(cat <<FEATURECOUNTS
echo "START FEATURECOUNTS"
cp \${workingdir}/3_GSnap/* \${localdir1}
cp \${refcountpath}/\${refcountfile} \${localdir1}

  echo "Subreading feature counts..."
   time featureCounts\
 	-a \${refcountfile}\
 	-C\
 	-Q 0\
 	-p\
 	-P\
 	-d 50\
 	-D 600\
 	-B\
 	-t \${gtftype}\
 	-F GTF\
 	-g \${gtfname}\
 	-o \${sortedbam}.counts.txt\
 	-T 15\
 	-s 2\
 	\${sortedbam}

## Copy output files back to the user\'s specified destination
echo "Copying output files back to your specified destination..."
time cp \${sortedbam}.counts.txt \${featureresults}
FEATURECOUNTS
)

COMBINE_FEATURES_4_1=$(cat <<COMBINE_FEATURES
echo "START COMBINING BAM FEATURECOUNTS"
cp \${workingdir}/3_GSnap/* \${localdir1}
cp \${refcountpath}/\${refcountfile} \${localdir1}
 echo "Combining subreading feature counts..."
   time featureCounts\
 	-a \${refcountfile}\
 	-C\
 	-Q 0\
 	-p\
 	-P\
 	-d 50\
 	-D 600\
 	-B\
 	-t \${gtftype}\
 	-F GTF\
 	-g \${gtfname}\
 	-o ALL.counts.txt\
 	-T 15\
 	-s 2\
 	*.sorted.bam

echo "Copying output files back to your specified destination..."
time cp ALL.counts.txt \${featureresults}
COMBINE_FEATURES
)

function GET_STEPS {
  if [ "${STEPS}" == "all" ]
  then
    echo "${FASTQC_1}"
    echo "${TRIM_2}"
    echo "${GSNAP_3}"
    echo "${FEATURECOUNTS_4}"
  elif [ "${STEPS}" == "1" ]
  then
    echo "${FASTQC_1}"
  elif [ "${STEPS}" == "2" ]
  then
    echo "${TRIM_2}"
  elif [ "${STEPS}" == "3" ]
  then
    echo "${GSNAP_3}"
  elif [ "${STEPS}" == "4" ]
  then
    echo "${FEATURECOUNTS_4}"
  elif [ "${STEPS}" == "12" ]
  then
    echo "${FASTQC_1}"
    echo "${TRIM_2}"
  elif [ "${STEPS}" == "13" ]
  then
    echo "${FASTQC_1}"
    echo "${GSNAP_3}"
  elif [ "${STEPS}" == "14" ]
  then
    echo "${FASTQC_1}"
    echo "${FEATURECOUNTS_4}"
  elif [ "${STEPS}" == "23" ]
  then
    echo "${TRIM_2}"
    echo "${GSNAP_3}"
  elif [ "${STEPS}" == "24" ]
  then
    echo "${TRIM_2}"
    echo "${FEATURECOUNTS_4}"
  elif [ "${STEPS}" == "34" ]
  then
    echo "${GSNAP_3}"
    echo "${FEATURECOUNTS_4}"
  elif [ "${STEPS}" == "4_1" ]	#Used for a separate run through
  then
  	echo "${COMBINE_FEATURES_4_1}"
  fi
}

#function EMAIL_TO {
#echo "EMAILING REPORT"

#if [ grep ${workingdir\${filename}-Result.o* ]

#}

for (( i=1; i<=$#; i++ ))
do
## Evaluates and registers the command line arguments
  source $(eval echo \${$i})

if [[ "${WORKING_DIR}" == "." ]]
then
  WORKING_DIR="$(pwd)"
else
  cd "${WORKING_DIR}"
fi

localdir1=/scratch3/$USER

## Determining file names
DAT_DIR=${SOURCEPATH}/${SOURCEDIR}/raw_data

ls -1 $DAT_DIR | sed -e 's/_R1.*//' -e 's/_R2.*//' | sort -u > filenames.txt

## Going to the location of the scripts - prints Palmetto Outputs here
#mv filenames.txt ${SOURCEPATH}/${SOURCEDIR}/${LOC_SCRIPT}
#mv filenames.txt ${localdir1}

## Creating Subfolders for Each Step
mkdir 1_FastQC
mkdir 2_Trim
mkdir 2.5_PostTrimFastQC
mkdir 3_GSnap
mkdir 4_FeatureCounts
mkdir Palmetto_Outputs

  ## Loop to create x amount of jobs
  tail filenames.txt | while read FILENAME
  do
    
    ##////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ## Starting job script
    ##////////////////////////////////////////////////////////////////////////////////////////////////////////////
    count=1
    
    JOB_SCRIPT="${FILENAME}-RNAseq${count}"
    (
    echo "#!/bin/bash"
    echo "#PBS -N ${FILENAME}-Result"
    echo "#PBS -l select=1:ncpus=8:mem=70gb:interconnect=10ge,walltime=32:00:00"
    echo "#PBS -j oe"
    echo "#PBS -m abe"
    echo "#PBS -q workq"
    echo ""
    #echo "set -x"
    echo "qstat -xf \${PBS_JOBID}"
    echo ""
    echo "${MODULES_NEEDED}"
    echo ""
    echo "workingdir=\"${WORKING_DIR}\" # absolute path"
    echo "datadir=\"${DAT_DIR}\" # absolute path"
    echo "fastqcdir=\"\${workingdir}/1_FastQC\""
    echo "trimsource=\"\${workingdir}/2_Trim\""
    echo "posttrimdir=\"\${workingdir}/2.5_PostTrimFastQC\""
    echo "resultsdir=\"\${workingdir}/3_GSnap\""
    echo "featureresults=\"\${workingdir}/4_FeatureCounts\""
    echo "filename=\"${FILENAME}\""
    echo "filename_forward=\"\${filename}_R1.paired.fastq.gz\""
    echo "filename_reverse=\"\${filename}_R2.paired.fastq.gz\""
    echo "sortedbam=\"\${filename}.sorted.bam\""
    echo "refname=\"${REFNAME}\""
    echo "refcountfile=\"${REFCOUNTFILE}\""
    echo "gtftype=\"${GTFTYPE}\""
    echo "gtfname=\"${GTFNAME}\""
    echo "loc_script=\"${LOC_SCRIPT}\""
    echo "refcountpath=\"${REFCOUNTPATH}\""
    echo "email=\"${EMAIL}\""
    echo ""
    echo "localdir1=\"/scratch3/$USER\""
    echo ""
    echo "mkdir -p \"\${localdir1}\""
    echo "cd \${localdir1}"
    echo "pwd"
    echo ""
    echo "## Copying raw data files to the local_scratch"
    echo "echo \"Copying raw data files to the local_scratch...\""
    echo "cp \${datadir}/\${filename}*.fastq.gz \${localdir1}"
    echo ""
    GET_STEPS
    echo ""
    echo "## Removing all unnecessary files"
    echo "echo \"Removing all unnecessary files..\""
    echo "rm \${resultsdir}/\${filename}*unpaired.fastq.gz"
    echo "rm \${resultsdir}/\${filename}*fastqc.zip"
    echo "rm \${resultsdir}/\${filename}*fastq.gz"
    echo "rm \${resultsdir}/\${filename}*.html"
    echo "rm \${featureresults}/\${filename}*fastqc.zip"
    echo "rm \${featureresults}/\${filename}*unpaired.fastq.gz"
    echo "rm \${featureresults}/\${filename}*fastq.gz"
    echo "rm \${featureresults}/\${filename}*.sam"
    echo "rm \${datadir}/\${filename}*.html"
    echo "rm \${datadir}/\${filename}*.zip"
    echo "rm \${trimsource}/\${filename}*.html"
    echo "rm \${trimsource}/\${filename}*.zip"
    echo "rm \${featureresults}/trimmomatic-0.36.jar"
    echo "rm \${featureresults}/illuminaClipping.fa"
    echo "rm \${featureresults}/\${filename}*.bai"
    echo "rm \${featureresults}/\${filename}*.html"
    echo "rm \${posttrimdir}/\${filename}*_00*"
    echo ""
    echo "echo \"Moving palmetto output files...\""
    #EMAIL_TO
    echo "mv \${workingdir}/*.o* \${workingdir}/Palmetto_Outputs"
    echo "mv \${workingdir}/*.pbs \${workingdir}/\${loc_script}"
    echo ""
    echo "## Changing Permissions"
    echo "chmod -R 770 \${workingdir}/*"
    echo ""
    echo "qstat -xf \"\${PBS_JOBID}\""
    echo ""
    ) > ${JOB_SCRIPT}.pbs
    echo "Created Host Watcher Script"
    
    ((count++))
    # currently in WORKING_DIR so all PBS output files will be here
    qsub ${JOB_SCRIPT}.pbs
    
    ##////////////////////////////////////////////////////////////////////////////////////////////////////////////
  done
done

echo "All jobs created"
