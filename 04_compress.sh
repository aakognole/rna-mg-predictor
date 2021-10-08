#!/bin/bash

for i in $(seq 1 1 5)
do
    echo ${i}
    tar -czf ${i}.tgz ./${i}
    rm -rf ${i}
done
