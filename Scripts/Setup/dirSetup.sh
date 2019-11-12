#!/bin/bash

if [[ ! -d $1 -a ! -h $1 ]]; then
    mkdir -p $1 || \
        { echo "$1 exists and is not a symlink or directory!" && false }
fi
