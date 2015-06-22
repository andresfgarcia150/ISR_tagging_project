#! /bin/bash

mensaje='Changing the folder names'
echo "Hi, $mensaje"

# loop over 100K event folders
head_folder='_Tops_result_WI'
sim_folder='_Tops_MG_1K_AG_'
max_val=7
for ((i=0; i<=${max_val}; i++))
do
        for ((j=($i+1); j<=${max_val}; j++))
        do
                for ((k=($j+1); k<=${max_val}; k++))
                do
                         mv ${head_folder}/Overall_sTopsWISR_${i}_${j}_${k}.txt ${head_folder}/_Tops_WI_Overall_${i}_${j}_${k}.txt
                done
        done
done
