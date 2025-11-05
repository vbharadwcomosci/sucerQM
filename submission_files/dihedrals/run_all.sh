##This script submits the dihedral scan calculations
#!/bin/sh
for parent in */; do
cd $parent
for dir in */; do
cd $dir
for i in */; do
cd $i
rm *.out* *clean* *xyz
cp ../../../gaussian_submission_kestrel.sh .
sbatch gaussian_submission_kestrel.sh
cd ../
done
cd ../
done
cd ../
done
