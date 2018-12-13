##############################################
# By: Jessica Qiu
# Last Update: December 13, 2017
# Purpose: Format Results
# How to Use: ./GO_term_compressor.sh <file>
##############################################

#!/bin/bash

awk -v FORMAT=," " '!tmp[$1,$2]++{arr[$1] =($1 in arr ? arr[$1] FORMAT : "" ) $2}END{for(i in arr)print i" "arr[i]}' $1 > format_output.csv
