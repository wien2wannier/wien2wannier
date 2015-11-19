#!/bin/bash
#stores directory name to filename variable 10/01/10

filename=$(basename $(pwd))
echo "seedname is: $filename"
export filename
