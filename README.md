# SUCER-qm 
### Saccharide pUcker Conformation and Energy exploreR using quantum mechanics

## Overview
- sucer-QM is a set of python codes that enable the generation and energetic exploration of puckered conformations for any 5-membered (furanose) or 6-membered (pyranose) monosaccharide structures.
- It accepts smile strings as inputs and generates puckered conformation files in pdb/xyz or gaussian input format based on specified cremer-pople parameters (q, $\theta$ for furanoses and q, $\theta$ and $\phi$ for pyranoses)
- Ensuring energy minima for each pucker state involves the full exploration of the orientation of all exocyclic groups attached to the saccharide. To enable this sucer-QM has pre-defined functions (inital_dihedral_scan) that enable QM energetic scans of dihedral angles of exocyclic groups. This is undertaken in an iterative manner with (repeat_dihedral_scan) functions, wherein each exocylcic group is scanned independently until no further lower-minima are observed.
- Finally, the final_unconstrained  or final_constrained functions enable the calculation of final optimized structures for thermodynamic and energetic calculations.

 Contributor: [Sean Brooks](https://github.com/sbrookse)

## Dependencies

This code requires the following pre-requisite packages to be installed.
- RDKit
- numpy
- pandas
- Matplotlib
- RingAnalysis (https://github.com/lucianlschan/RING)
- RingReconstruction (https://github.com/lucianlschan/RING)
- py_rdl (RingDecomposerLib) (https://github.com/rareylab/RingDecomposerLib)


## Running sucer-QM
* Step 0
  - Update the smile strings of the sugar molecules whose pucker states are desired to be explored
  - Edit the pucker-parameters.csv file(s) to ensure that the parameters (esp. q values) correspond to the desired pucker states to be explored.
* Step 1 - Run initial optimizations for each sugar and its various pucker states
  - ```
    python
    import pucker_gen as pg
    pg.generate_puckers('input')
    ```
    - This will generate a folder called input containing folders named after each sugar whose pucker state is desired to be explored.
    - Each sugar folder contains sub-folders named after the pucker state and contains gaussian input files and .xyz files of the pucker state
  - If you have access to a HPC workstation with gaussian installed then you may run the following command to submit the QM calcualtions
  - ```
    sh submit.sh
    sh compile.sh
    ```
    - 
* Step 2 - Perform dihedral scans for exocyclic groups attached to the sugar
  - kdjfa'lsdf
  - ```
    ksjdfasdf
    asdfasdf
    ```
    - This will generate gaussian input files for running dihedral scans for each individual exocyclic group.
    -
* Step 3 - Perform final optimizations of the lowest energy configurations identified from Step 2
  - asdfkjasdf
  - asdfkljasldkfj
  - 
aaaa
