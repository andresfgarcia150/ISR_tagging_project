#! /bin/bash

mensaje='Executing several matching procedures'
echo "Hi, $mensaje"

# loop over 100K event folders
exe_file='ISR_matching'
conf_file='config_file.txt'
for i in {261..270}
do
	./$exe_file $conf_file $i
done

