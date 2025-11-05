#!/usr/bin/env python                                                                                                 
# python Extract_Optimized_From_Gaussian.py filename                                                                    
import sys, os
from pathlib import Path

def extract_all(text):
    opt_param_line_nums = []
    input_orient_line_nums = []
    opt_scf = []
    end_coor_line_nums = []
    scf_temp = 0
    # Start count at 1 because files start at line 1 not 0                                                              
    count = 1
    if "Normal termination" in text[len(text)-1]:
        for line in text:
            if 'SCF Done' in line:
                scf_temp = line.split()[4]
            elif 'Optimized Parameters' in line:
                opt_param_line_nums.append(count)
                opt_scf.append(scf_temp)
            elif 'Input orientation' in line:
                input_orient_line_nums.append(count)
            elif '---------' in line:
                end_coor_line_nums.append(count)
            count += 1
    else:
        opt_scf.append('N/A')

    # Add 5 because the first actual coordinate is +5 lines                                                            
    input_orient_line_nums = [x+5 for x in input_orient_line_nums]
    
    return opt_param_line_nums, opt_scf, input_orient_line_nums, end_coor_line_nums

def write_xyz(text, output_folder, sugar, angle, substituent, count):
    code = {"1" : "H", "2" : "He", "3" : "Li", "4" : "Be", "5" : "B", \
            "6"  : "C", "7"  : "N", "8"  : "O", "9" : "F", "10" : "Ne", \
            "11" : "Na" , "12" : "Mg" , "13" : "Al" , "14" : "Si" , "15" : "P", \
            "16" : "S"  , "17" : "Cl" , "18" : "Ar" , "19" : "K"  , "20" : "Ca", \
            "21" : "Sc" , "22" : "Ti" , "23" : "V"  , "24" : "Cr" , "25" : "Mn", \
            "26" : "Fe" , "27" : "Co" , "28" : "Ni" , "29" : "Cu" , "30" : "Zn", \
            "31" : "Ga" , "32" : "Ge" , "33" : "As" , "34" : "Se" , "35" : "Br", \
            "36" : "Kr" , "37" : "Rb" , "38" : "Sr" , "39" : "Y"  , "40" : "Zr", \
            "41" : "Nb" , "42" : "Mo" , "43" : "Tc" , "44" : "Ru" , "45" : "Rh", \
            "46" : "Pd" , "47" : "Ag" , "48" : "Cd" , "49" : "In" , "50" : "Sn", \
            "51" : "Sb" , "52" : "Te" , "53" : "I"  , "54" : "Xe" , "55" : "Cs", \
            "56" : "Ba" , "57" : "La" , "58" : "Ce" , "59" : "Pr" , "60" : "Nd", \
            "61" : "Pm" , "62" : "Sm" , "63" : "Eu" , "64" : "Gd" , "65" : "Tb", \
            "66" : "Dy" , "67" : "Ho" , "68" : "Er" , "69" : "Tm" , "70" : "Yb", \
            "71" : "Lu" , "72" : "Hf" , "73" : "Ta" , "74" : "W"  , "75" : "Re", \
            "76" : "Os" , "77" : "Ir" , "78" : "Pt" , "79" : "Au" , "80" : "Hg", \
            "81" : "Tl" , "82" : "Pb" , "83" : "Bi" , "84" : "Po" , "85" : "At", \
            "86" : "Rn" , "87" : "Fr" , "88" : "Ra" , "89" : "Ac" , "90" : "Th", \
            "91" : "Pa" , "92" : "U"  , "93" : "Np" , "94" : "Pu" , "95" : "Am", \
            "96" : "Cm" , "97" : "Bk" , "98" : "Cf" , "99" : "Es" ,"100" : "Fm", \
            "101": "Md" ,"102" : "No" ,"103" : "Lr" ,"104" : "Rf" ,"105" : "Db", \
            "106": "Sg" ,"107" : "Bh" ,"108" : "Hs" ,"109" : "Mt" ,"110" : "Ds", \
            "111": "Rg" ,"112" : "Uub","113" : "Uut","114" : "Uuq","115" : "Uup", \
            "116": "Uuh","117" : "Uus","118" : "Uuo"}

    opt_param_line_nums, opt_scf, input_orient_line_nums, end_coor_line_nums = extract_all(text)
    energy_out = open(output_folder+'/outputs/all_energies.txt', 'a')
    
    energy_out.write(f'Sugar: {sugar}\n')
    energy_out.write('----------------------\n')
    energy_out.write('----------------------\n')

    if count == 1:
        energy_out.write('Pucker State: '+pucker+'\n')
        energy_out.write('----------------------\n')
        
    if opt_scf[0] != 'N/A':
        new_input_lns = []
        count = 0
        for x in reversed(input_orient_line_nums):
            if count > len(opt_param_line_nums)-1:
                break
            if x < list(reversed(opt_param_line_nums))[count]:
                    new_input_lns.append(x)
                    count += 1

        new_input_lns.reverse()

        new_end_lns = []
        count = 0
        for x in end_coor_line_nums:
            if count > len(new_input_lns)-1:
                break
            if x > new_input_lns[count]:
                    new_end_lns.append(x)
                    count += 1

        intervals = zip(new_input_lns, new_end_lns)

        # Convert each interval into a .xyz                                                                     

        for i, intvl in enumerate(intervals):
            n_atoms = intvl[1] - intvl[0]
            angle = i * 15
            xyz_name = output_folder+'/outputs/xyz_outputs/'+pucker+'_'+substituent+'_'+str(angle)+'.xyz'
            #pdb_name = output_folder+'/outputs/pdb_outputs/'+pucker+'_'+substituent+'_'+str(angle)+'.pdb'
            xyz_out = open(xyz_name, 'w')
            xyz_out.write(str(n_atoms)+'\n\n')
            for x in range(intvl[0]-1, intvl[1]-1):
                column = text[x].split()
                print(code[column[1]],float(column[3]),float(column[4]),float(column[5]),file=xyz_out)
            xyz_out.close()
            #os.system(f'obabel -ixyz {xyz_name} -opdb -O {pdb_name}')

            energy_out.write('Substituent: '+substituent+'\n')
            energy_out.write('Angle: '+str(angle)+'\n')
            energy_out.write('Energy: '+str(opt_scf[i])+'\n\n\n')
    else:
        energy_out.write('Substituent: '+substituent+'\n')
        energy_out.write('Angle: 0\n')
        energy_out.write('Energy: N/A\n\n\n')
    energy_out.close()

root = Path.cwd()
root_path = Path(root)

for sugar in os.listdir('.'):
    sugar_path = os.path.join(root, sugar)
    if os.path.isdir(sugar_path):
        os.system(f'rm -r ./{sugar}/outputs')
        os.system(f'mkdir ./{sugar}/outputs')
        os.system(f'mkdir ./{sugar}/outputs/xyz_outputs')
        #os.system(f'mkdir ./{sugar}/outputs/pdb_outputs')
        for pucker in os.listdir(sugar_path):
            pucker_path = os.path.join(sugar_path, pucker)
            if os.path.isdir(pucker_path) and not 'outputs' in pucker:
                count = 0
                substituents = os.listdir(pucker_path)
                sorted_substituents = sorted(substituents, key=lambda x: x)
                for substituent in sorted_substituents:
                    substituent_path = os.path.join(pucker_path, substituent)
                    if os.path.isdir(substituent_path):
                        count += 1
                        filename = pucker+'.log'
                        logfile = os.path.join(substituent_path, filename)
                        logfile_bn = os.path.splitext(logfile)[0]
                        logfile_fh = open(logfile, 'r')
                        text = logfile_fh.readlines()
                        logfile_fh.close()
                        write_xyz(text, sugar_path, sugar, pucker, substituent, count)
