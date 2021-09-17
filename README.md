# Molecule-Adsorption

This repository includes 3 example folders to employ the .py scripts to generate a .txt file with the main features describing the adsorption modes of molecules in atomic clusters.

TMx_CH3_or_CH4: This folder includes the python script ch43_new2.py, which can be runned simply with "python3 ch43_new2.py". It will consider all the .xyz files in the same folder and generate a table (table_feat.txt) with the main adsorption features. Its work for any dataset of metallic clusters (TMx, x = 1, 2, 3, ... n.) interacting with a CH4 or a CH3 - This program identifies automatically weather it is a CH4 or CH3 interacting with the atomic cluster by counting the number of hydrogens in the .xyz file, be careful.

TMx_CH3+H: This folder includes the python script ch3h_new2.py, which can be runned simply with "python3 ch3h_new2.py". It will consider all the .xyz files in the same folder and generate a table (table_feat.txt) with the main adsorption features. Its work for any dataset of metallic clusters (TMx, x = 1, 2, 3, ... n.) interacting with a CH3+H.

TMx_H: This folder includes the python script htm.py, which can be runned simply with "python3 htm.py". It will consider all the .xyz files in the same folder and generate a table (table_feat.txt) with the main adsorption features. Its work for any dataset of metallic clusters (TMx, x = 1, 2, 3, ... n.) interacting with an H. It considers that only one hydrogen is present in the .xyz file, thus, be careful.
