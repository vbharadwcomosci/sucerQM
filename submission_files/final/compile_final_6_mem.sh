##This code does the following
## - creates an output folder and organize all outputs from the final optimization                       
## - generates pdb and xyz () files from gaussian log files for the next set of calcs                      
## - outputs all_energies.txt containing all energies, and pucker parameters for the optimized geoms       
##   calculated using cpptraj
#!/bin/sh
module load amber
${pucker}_opt.xyz
for parent in */; do
cd $parent
sugar="${parent::-1}"
rm -rf outputs
mkdir outputs
touch outputs/all_energies.txt
mkdir outputs/pdb_outputs
touch outputs/pdb_outputs/opt_all.pdb
touch outputs/pdb_outputs/ptraj.in
mkdir outputs/xyz_outputs
echo "trajin outputs/pdb_outputs/opt_all.pdb" > outputs/pdb_outputs/ptraj.in
echo "pucker outputs/pdb_outputs/opt_all @2 @3 @4 @5 @6 @1 range360 theta amplitude cremer out outputs/pdb_outputs/pucker_check.dat" >>\
 outputs/pdb_outputs/ptraj.in
echo "Ring Atoms: 6" >> outputs/all_energies.txt
echo "Sugar: $sugar" >> outputs/all_energies.txt
echo -e "--------------------------------\n--------------------------------" >> outputs/all_energies.txt

for i in */; do
if [[ ! "$i" == *"output"* && ! "$i" == *"input"* ]]; then
cd $i
pucker="${i::-1}"
echo "Sugar: $sugar -- Pucker State: $pucker"
if tac *.log | grep -q "Normal termination"; then
obabel ${pucker}.log -oxyz -O ${pucker}_opt.xyz
obabel -ixyz ${pucker}_opt.xyz -opdb > ../outputs/pdb_outputs/${pucker}_opt.pdb
cat ../outputs/pdb_outputs/${pucker}_opt.pdb >> ../outputs/pdb_outputs/opt_all.pdb
cp ${pucker}_opt.xyz ../outputs/xyz_outputs
electronic=$(grep "Energy" ${pucker}_opt.xyz)
else
echo "!!!!RUN FAILED!!!!"
fi
cd ../
fi
done

cpptraj -i outputs/pdb_outputs/ptraj.in -p outputs/pdb_outputs/14B_opt.pdb

i=1
for dir in */; do
pucker="${dir::-1}"
if [[ ! "$dir" == *"output"* && ! "$dir" == *"input"* ]]; then
cd $dir
echo "Pucker State: $pucker" >> ../outputs/all_energies.txt
echo "--------------------------------" >> ../outputs/all_energies.txt
if tac *.log | grep -q "Normal termination"; then
((i++))
ZPE_correction_line=$(grep "Zero-point correction" *.log)
electronic_ZPE_line=$(grep "Sum of electronic and zero-point Energies" *.log)
electronic_thermal_energy_line=$(grep "Sum of electronic and thermal Energies" *.log)
electronic_thermal_enthalpy_line=$(grep "Sum of electronic and thermal Enthalpies" *.log)
electronic_thermal_FE_line=$(grep "Sum of electronic and thermal Free Energies" *.log)

ZPE_correction=$(echo "$ZPE_correction_line" | awk '{print $3}')
electronic_ZPE=$(echo "$electronic_ZPE_line" | awk '{print $7}')
electronic_energy=$(echo "$ZPE_correction + $electronic_ZPE" | bc)
electronic_thermal_energy=$(echo "$electronic_thermal_energy_line" | awk '{print $7}')
electronic_thermal_enthalpy=$(echo "$electronic_thermal_enthalpy_line" | awk '{print $7}')
electronic_thermal_FE=$(echo "$electronic_thermal_FE_line" | awk '{print $8}')

pucker_line=$(sed -n "${i}p" "../outputs/pdb_outputs/pucker_check.dat")
phi=$(echo "$pucker_line" | awk '{print $2}')
q=$(echo "$pucker_line" | awk '{print $3}')
theta=$(echo "$pucker_line" | awk '{print $4}')
else
electronic="N/A"
phi="N/A"
q="N/A"
theta="N/A"
fi
echo "Phi: $phi" >> ../outputs/all_energies.txt
echo "Q: $q" >> ../outputs/all_energies.txt
echo "Theta: $theta" >> ../outputs/all_energies.txt
echo -e "Electronic Energy: $electronic_energy" >> ../outputs/all_energies.txt
echo -e "Sum of electronic and zero-point Energies: $electronic_ZPE" >> ../outputs/all_energies.txt
echo -e "Sum of electronic and thermal Energies: $electronic_thermal_energy" >> ../outputs/all_energies.txt
echo -e "Sum of electronic and thermal Enthalpies: $electronic_thermal_enthalpy" >> ../outputs/all_energies.txt
echo -e "Sum of electronic and thermal Free Energies: $electronic_thermal_FE\n\n" >> ../outputs/all_energies.txt
cd ../
fi
done
cd ../
done
