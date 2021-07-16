#!/bin/bash

# 2021.7.16 TLCL 16

# Check the network settings
#$ netstat [options]
## options
#-ie network interface
#-r display the kernel’s network routing table

netstat -ie
netstat -r
