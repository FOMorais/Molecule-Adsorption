# Molecule-Adsorption

This repository includes 3 example folders to employ the .py scripts to generate a .txt file with the main features describing the adsorption modes of molecules in atomic clusters.

--TMx_CH3_or_CH4: This folder includes the python script ch43_new2.py, which can be runned simply with "python3 ch43_new2.py". It will consider all the .xyz files in the same folder and generate a table (table_feat.txt) with the main adsorption features. Its work for any dataset of metallic clusters (TMx, x = 1, 2, 3, ... n.) interacting with a CH4 or a CH3 - This program identifies automatically weather it is a CH4 or CH3 interacting with the atomic cluster by counting the number of hydrogens in the .xyz file, be careful.

The features calculated for TMxCH3 or TMxCH4 follow this order: filename, average HCH angle, distance of the closest H (Hc) to the closest TM (TMc) to the molecule, distance of the farthest H (Hf) to the closest TM (TMc) to the molecule, distance of the C to the closest TM (TMc) to the molecule, dihedral angle between Hf/C/Hc/TMc, coordination of the closest TM (TMc) to the molecule, average bond length of the closest TM (TMc) to the molecule, orientation angle between Hc/C/TMc, coordination of the C.

--TMx_CH3+H: This folder includes the python script ch3h_new2.py, which can be runned simply with "python3 ch3h_new2.py". It will consider all the .xyz files in the same folder and generate a table (table_feat.txt) with the main adsorption features. Its work for any dataset of metallic clusters (TMx, x = 1, 2, 3, ... n.) interacting with a CH3+H.

The features calculated for TMxCH3+H follow this order: filename, average HCH angle, distance of the closest H (Hc) to the closest TM (TMc) to the molecule, distance of the farthest H (Hf) to the closest TM (TMc) to the molecule, distance of the C to the closest TM (TMc) to the molecule, dihedral angle between Hf/C/Hc/TMc, coordination of the closest TM (TMc) to the molecule, average bond length of the closest TM (TMc) to the molecule, orientation angle between Hc/C/TMc, coordination of the C.

file1.write(str(filename)+"\t"+str(round(HCHav,5))+"\t"+str(round(DHcTM,5))+"\t"+str(round(DHfTM,5))+"\t"+str(round(DCtm,5))+"\t"+str(round(abs(CHads),5))+"\t"+str(round(abs(Hdied),5))+"\t"+str(round(abs(HCTM),5))+"\t"+str(round(abs(CNTM),5))+"\t"+str(round(abs(davTM),5))+"\t"+str(round(abs(CNHads),5))+"\t"+str(round(abs(davHads),5))+"\t"+str(round(abs(CNC),5))+"\n")

TMx_H: This folder includes the python script htm.py, which can be runned simply with "python3 htm.py". It will consider all the .xyz files in the same folder and generate a table (table_feat.txt) with the main adsorption features. Its work for any dataset of metallic clusters (TMx, x = 1, 2, 3, ... n.) interacting with an H. It considers that only one hydrogen is present in the .xyz file, thus, be careful.
