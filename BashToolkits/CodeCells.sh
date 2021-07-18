#!/bin/bash


#$ echo [string]

echo test is good  # even without the quotes
# echo usage
echo 'Mesages for users'now=$("date")
echo -e 'Date: \n $now'
echo -e "Date: \n $now"
echo -e "Date: \n \$now"

## Options
#-e  # enable interpretation of following backslash escapes