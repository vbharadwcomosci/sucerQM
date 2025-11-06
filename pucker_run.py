##This code will run pucker_gen.py to perform the sugar pucker conformatoinal and energetic exploration
##In Step 1 we generate the puckered conformations for the desired saccharide
import os as os
import pucker_gen as pg
pg.generate_puckers('step1')

##This will generate a folder step1 and place all the files in there
##At this stage one may submit the calculations on a cluster using the run_all.sh script in the 5_mem/6_mem folder created in step1
##After the calculations finish, run the compile.sh script in the 5_mem/6_mem folder

##In step 2 we run dihedral scans for each pucker state to ensure we are at a 'global' minima for that pucker state
##This will generate a folder called step2 in the current directory.
##This will again contain step2/5mem(or 6mem)/sugar/puckerstate/dihedral

#pg.inital_dihedral_scan('step2','step1')

##Now one can run the dihedral scans run_all.sh script as in step1

##To fix any convergence issues, one can use the ./rerun_failed.sh or ./rerun_failed2.sh or ./rerun_failes3.sh to resubmit scans.

##Step3...N
##In the case of any new minima found from step2, one should re-run the other dihedrals to ensure no lower minima are accessible.
#pg.repeat_dihedral_scan('step3','step2')

##Perform as many repeat scans as required.

##StepN+1
##Once a 'global' minima for each sugar and pucker state is obtained, one can run constrained optimizations at desired level of theory 
##and basis-set and consider solvations etc. To do this....
## a. prepare a folder for each sugar
## b. make an all inputs folder and copy all the finished scans from /step3(N)/5mem(or 6mem)/sugar/outputs/Finished_Structures/
##Needs to be run for each sugar individually

#pg.final_constrained('stepN+1',5(or 6))

##StepN+2
#pg.final_unconstrained('stepN+2')
