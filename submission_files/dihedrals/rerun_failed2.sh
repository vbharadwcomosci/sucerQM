## This code is to be run after initial dihedral scans complete
## checks for normal termination
## if any of the optimizations fail then this code edits the dihedral scan input files to go the
## original direction but changes the dihedral constraint on the ring and resubmits the scans
##
#!/bin/sh
for parent in */; do
cd $parent
for dir in */; do
if [[ ! "$dir" == *"output"* && ! "$dir" == *"input"* ]]; then
cd $dir
for i in */; do
cd $i
if ! tail -n 1 *.log | grep -q "Normal termination"; then
sed -i 's/D 2 3 4 5 F/D 3 4 5 6 F/' *.com
sed -i 's/S 24 -15.0/S 24 15.0/' *.com
rm *.out* *clean* *xyz
cp ../../../gaussian_submission_kestrel.sh .
sbatch gaussian_submission_kestrel.sh
fi
cd ../
done
cd ../
fi
done
cd ../
done
