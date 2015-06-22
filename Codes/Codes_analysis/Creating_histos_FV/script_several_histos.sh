#! /bin/bash

mensaje='Executing several matching procedures'
echo "Hi, $mensaje"

# loop over 100K event folders
exe_file='Creating_histo'
max_val=7
for ((i=0; i<=${max_val}; i++))
do
	for ((j=($i+1); j<=${max_val}; j++))
	do
		for ((k=($j+1); k<=${max_val}; k++))
		do
			./$exe_file $i $j $k
		done
	done
done

