#!/bin/bash

if [ ! -e toolset ]; then
    wget https://github.com/mmtsb/toolset/archive/refs/heads/master.zip
    unzip -q master.zip
    rm master.zip
    mv toolset-master toolset
else
    exit
fi
