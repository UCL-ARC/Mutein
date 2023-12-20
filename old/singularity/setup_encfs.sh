#!/bin/bash

#check correct working directory
mutein check-cwd || { echo "check-cwd failed"; exit 1; }

#create encfs key in encrypted file
mutein encrypt-key --output-file config/encfs_password || { echo "encrypt-key failed"; exit 1; }

#ensure cipher and plain directories exist
mkdir -p cipher plain || { echo "creating folders failed"; exit 1; }
chmod go-wrx cipher plain || { echo "setting permissions failed"; exit 1; }