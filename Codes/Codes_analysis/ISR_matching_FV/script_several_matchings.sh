#! /bin/bash

mensaje='Executing several matching procedures'
echo "Hi, $mensaje"

# loop over 100K event folders
exe_file='ISR_matching'
for i in {001..010}
do
	./$exe_file $i
done

