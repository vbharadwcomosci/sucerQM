#!/bin/sh
for dir in */; do
cd $dir
for i in */; do
cd $i
rm *.log* *chk* *out *err
cp ../../gaussian_submission_kestrel.sh .
sbatch gaussian_submission_kestrel.sh
cd ../
done
cd ../
done
