#! /bin/bash

mensaje='Executing several matching procedures'
echo "Hi, $mensaje"

# loop over 100K event folders
exe_file='ISR_matching_FV'
for i in {000..009}
do
	./$exe_file $i
done

