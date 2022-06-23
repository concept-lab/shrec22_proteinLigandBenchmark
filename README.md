# shrec22_proteinLigandBenchmark
Dataset and evalutaion tools of the Shrec 2022  contest on protein-ligand binding site recognition

## Usage 
For evaluating putative pocket in PQR or OFF format:

python3 evaluate.py \<directoryName containing participant results\>

### Extra

The script filterLig.py, creates a folder containing all ligands (in xyz format) keeping only ligand atoms within 5A from any correspondent protein atom.
These are actually the ligand coordinates used for evaluation (the same filtering process is performed within evaluate.py by providing the structures pqrs contained in the "allStructures" folder).

### NOTE
A lighter version (without the full database and PQR structure files)of the contest's participants evaluation tool is provided in https://github.com/concept-lab/shrec22_PLBinding_evaluationTools.git
