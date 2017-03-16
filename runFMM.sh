#!/bin/bash

FILE=output.txt

nodes=(1 2) 
#3 4 5 6 7)

for n in ${nodes[@]}; do
	./bbfmm -l 2 -g 1 -o $n >> $FILE
done

