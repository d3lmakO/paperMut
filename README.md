# paperMut
C++ and bash scripts for the paper: Accurate quantification of mutation rates in patient-derived colorectal cancer organoids via sequencing-based fluctuation assays.


## Description

This is the main folder where you can find the subfolders containing various parts of the program that generate the final lists of potential mutated bases.

Each subfolder contains its own README file, which provides more detailed information on the input and output for each part of the program.

To obtain the final list, the various subprograms should be executed in the following order:

1. parsing

2. common

3. mean coverage

4. modal coverage

Once the four subprograms have been executed, the list of potential mutations for CRC or yeast can be obtained using the subprogram

6. list mutations 

If you want the coordinates of the potential mutations, you can run directly

7. list coord mut

in the **signatures** dir.

You can obtain the list of the potential signatures using the exec

8. signatures mut

in the **signatures** dir.

**CRC_clonal_mr** and **YEAST_clonal_mr** are used to compute clonal mutation rate.  
To calculate this mutation rate, one needs an experimental parameter: the number of _generations_ (or number of cell divisions) during which the mutation accumulation experiment took place.

The **yeast/** folder contains bash scripts that generate mutation lists for yeast genomes. These scripts call some parts of the programs listed above,  
so the Makefile containing the subprograms must be compiled before analyzing the yeast genomes.

The **simulation/** folder, on the other hand, contains all the files related to the simulation, which can be compiled and run separately from the other components of the data analysis.  
A dedicated Makefile is included in this folder.

## Compilation

To compile a single subprogram run:
```bash
make X
```
with X a name written in the Makefile located in the main folder (this folder).

To compile all subprograms run:
```bash
make all
```
If you want to clean object files and executable files run:
```bash
make clean
```
or:
```bash
make cX
```

## Execution

To run a single subprogram follow this (example with parsing):
```bash
cd parsing/exec/
./parsing
```
Sometime the executable file requires additional arguments (the description is provided in the individual README files).  
All programs have their own executable in the ``exec`` subfolder.

## Dependencies

- g++20
- boost::iostreams
- gnu gsl library (only for the simulation)
- bwa 0.7.17-r1188 (only for the yeasts)
- samtools 1.20 (only for the yeasts)

## Citation

**Paper:**
> Pompei, S., Geroldi, A., Rivetti, P., Grassi, E., Vurchio, V., Tallarico, G., Corti, G., Tattini, L., Liti, G., Bertotti, A., Cosentino Lagomarsino, M. μSeq: Universal mutation rate quantification via deep sequencing of a single clonal expansion.

**Code Archive:**
> Geroldi, A. Code archive for "$\mu$Seq: Universal mutation rate quantification via deep sequencing of a single clonal expansion" (paperMut Version v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.21480295
