#!/bin/bash

# This should fail
hapmixmap test.conf

# If it didn't fail, return an error
if [ "$?" = "0" ]
then
	exit 1
else
	exit 0
fi
