#! /bin/bash

mensaje='Changing the folder names'
echo "Hi, $mensaje"

# loop over 100K event folders
folder='_Tops_matchs_WI'
file_move='ISR_jets_Tops_'
temp_str='WI_'
ext='.bn'
for i in {000..260}
do
	mv ${file_move}$i${ext} ${folder}/${file_move}${temp_str}$i${ext}
done

