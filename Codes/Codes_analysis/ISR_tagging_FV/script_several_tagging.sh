#! /bin/bash

mensaje='Executing several tagging procedures'
echo "Hi, $mensaje"

# loop over 100K event folders
exe_file='ISR_tagging'
conf_file='config_file.txt'
max_val=7
for ((i=0; i<=${max_val}; i++))
do
	for ((j=($i+1); j<=${max_val}; j++))
	do
		for ((k=($j+1); k<=${max_val}; k++))
		do
			./$exe_file $conf_file $i $j $k
		done
	done
done

