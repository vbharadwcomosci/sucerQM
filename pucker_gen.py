import numpy as np
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import rdGeometry as geom
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdMolTransforms as rdmt
from rdkit.Chem import Draw
from rdkit.Chem import rdDetermineBonds
import os, sys, inspect
import glob
import Ring_Analysis as RA
import Ring_Reconstruction as RR
import py_rdl
import openpyxl
import pandas as pd
import re
os.chdir(os.path.dirname(__file__))

# Get 5- and 6-membered pucker info that is used to provide input to the RR/RA package
# These files are found in the same folder as the package when it is downloaded and can be edited
c5_param_path = 'c5_pucker_params.csv'
params = pd.read_csv(c5_param_path)
c5_params = [params['q2'], params['Phi_radians']]
c5_states = params['State']
c6_param_path = 'c6_pucker_params.csv'
params = pd.read_csv(c6_param_path)
c6_params = [params['q2'], params['q3'], params['Phi_radians']]
c6_states = params['State']
sugar_info_path = 'Sugar_Info.xlsx'
sugar_info = pd.read_excel(sugar_info_path)

# Circularly traverse numbering list to put ring atoms in correct order starting with C1
def circular_traverse(ring_atoms, c1_index, o_index):
    n = len(ring_atoms)
    i = ring_atoms.index(c1_index)
    j = ring_atoms.index(o_index)

    # Determine direction for sort
    if ring_atoms[(i + 1) % n] == o_index:
        sorted_ring = [ring_atoms[(i - k) % n] for k in range(n)]
    else:
        sorted_ring = [ring_atoms[(i + k) % n] for k in range(n)]

    return sorted_ring

def renumber_sugar(SMILES, num_ring_atoms):
    mol = Chem.MolFromSmiles(SMILES)
    renumber = []
    ring_atoms = mol.GetRingInfo().AtomRings()[0]
    o_index = 0
    o_neighbors = ''
    all_subs = []

    # Find location of ring O and get neighboring C atoms
    for atom_index in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetSymbol() == 'O':
            o_index = atom.GetIdx()
            o_neighbors = atom.GetNeighbors()
            break

    # Find index of C1 to start numbering
    # C1 will be the one of the two ring O neighbors that has a neighboring non-ring O atom
    c1_index = 0
    for neighbor in o_neighbors:
        temp = neighbor.GetNeighbors()
        subs = [n for n in temp if n.GetIdx() not in ring_atoms]
        for atom in subs:
            if atom.GetSymbol() == 'O':
                c1_index = neighbor.GetIdx()

    circular = circular_traverse(ring_atoms, c1_index, o_index)
    if num_ring_atoms == 6:
        circular = [circular[-1]] + circular[:-1] 

    renumber = list(circular)
    dihedrals = []
    for atom_idx in circular:
        atom = mol.GetAtomWithIdx(atom_idx)

        for neighbor in atom.GetNeighbors():
            carboxylic = False
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    carboxylic = True

            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in circular:
                visited = []
                to_visit = [neighbor_idx]

                while to_visit:
                    curr = to_visit.pop()
                    visited.append(curr)

                    for neighbor2 in mol.GetAtomWithIdx(curr).GetNeighbors():
                        neighbor2_idx = neighbor2.GetIdx()
                        if neighbor2_idx not in circular and neighbor2_idx not in visited:
                            to_visit.append(neighbor2_idx)
                renumber.extend(visited)

                if len(visited) >= 2:
                    prev_atom_idx = circular[circular.index(atom_idx) - 1]
                    if not(carboxylic):
                        dihedrals.append([prev_atom_idx] + [atom_idx] + visited[0:2])
                    else:
                        dihedrals.append([prev_atom_idx] + [atom_idx] + [visited[0]] + [visited[2]])
                    if len(visited) > 2 and not(carboxylic):
                        dihedrals.append([atom_idx] + visited[0:3])


    mol_renumbered = Chem.RenumberAtoms(mol, renumber)
    new_numbering = []
    for atom in mol_renumbered.GetAtoms():
        new_numbering.append(atom.GetIdx())

    renumbered_dict = dict(zip(renumber, new_numbering))
    renumbered_dihedrals = []
    for dihedral in dihedrals:
        renumbered_dihedral = [renumbered_dict[dihedral_atom_idx] + 1 for dihedral_atom_idx in dihedral]
        renumbered_dihedrals.append(renumbered_dihedral)

    return mol_renumbered, renumbered_dihedrals

# Generates an rdkit Mol object from an xyz file
def get_mol_from_xyz(path):
    try:
        temp_mol = Chem.MolFromXYZFile(path)
        mol = Chem.Mol(temp_mol)
    except:
        with open(path, 'r') as file:
            xyz_string = file.read()

        pattern = r'\d+\.?\d*[eE][+-]?\d+'
        matches = re.finditer(pattern, xyz_string)
        replacements = []

        for match in matches:
            scientific_str = match.group(0)
            decimal_float = float(scientific_str)
            decimal_str = f"{decimal_float:.10f}".rstrip('0').rstrip('.')
            replacements.append((scientific_str, decimal_str))

        for old_val, new_val in replacements:
            xyz_string = xyz_string.replace(old_val, new_val)

        with open(path, 'w') as file:
            file.write(xyz_string)

        temp_mol = Chem.MolFromXYZFile(path)
        mol = Chem.Mol(temp_mol)

    rdDetermineBonds.DetermineBonds(mol)
    return mol

# Converts mol to xyz for writing to Gaussian input file
def get_xyz_from_mol(mol):
    xyz_output = Chem.MolToXYZBlock(mol)
    lines = xyz_output.splitlines()
    line_list = [line for i, line in enumerate(lines) if i > 1]
    output = "\n".join(line_list)
    return output

# Returns a specified dihedral angle
def get_dihedral(mol, dihedral_indices):
    conf = mol.GetConformer()
    dihedral_angle = rdMolTransforms.GetDihedralDeg(conf, *dihedral_indices)
    if dihedral_angle < 0:
        dihedral_angle += 360
    return dihedral_angle

# Returns a string of ring atom dihedral constraints for writing to Gaussian input file
def get_ring_dihedral_constraints(ring_atoms):
    if ring_atoms == 5:
        ring_dihedrals = 'D 1 2 3 4 F\nD 3 4 5 1 F'
    else:
        ring_dihedrals = 'D 1 2 3 4 F\nD 2 3 4 5 F\nD 4 5 6 1 F'
    return ring_dihedrals

# Returns a string of substituent dihedral atoms to constrain during the optimization
def get_scan_dihedral(dihedral_indices):
    scan_dihedral = 'D '
    for i in dihedral_indices:
        scan_dihedral = scan_dihedral + str(i) + ' '
    scan_dihedral += 'S 24 15.0'
    return scan_dihedral

# Optimizes only Hs on given molecule while keeping all other atoms in a fixed position
def optimize_H2(mol):
    mol_properties = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, mol_properties)

    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetSymbol() != 'H':
            ff.AddFixedPoint(i)

    ff.Minimize()
    return mol

# Finds lowest energy conformer out of 50 randomly generated conformers by minimizing RMSD
def align(base_mol, pucker_mol, ring_atoms):
    confs = Chem.Mol(pucker_mol)
    AllChem.EmbedMultipleConfs(confs, 50, randomSeed=0xf00d)

    atom_map = [[i,i] for i in range(ring_atoms)]
    rmsd = 100
    best_conf = 0

    for i in range(confs.GetNumConformers()):
        rmsd_temp = Chem.rdMolAlign.AlignMol(confs, base_mol, prbCid=i, atomMap=atom_map)
        if rmsd_temp < rmsd:
            rmsd = rmsd_temp
            best_conf = i

    for i in range(confs.GetNumConformers()):
        if i != best_conf:
            confs.RemoveConformer(i)
    return confs
        
# Writes geometry optimization input file for the initial structure
def write_gaussian_initial(mol, ring_atoms, state, output_path):
    xyz_output = get_xyz_from_mol(mol)
    ring_dihedrals = get_ring_dihedral_constraints(ring_atoms)
    
    with open(f'{output_path}/{state}/{state}.com', 'w') as f:
        f.write(f'''%chk={state}.chk
%nprocs=32
%mem=30GB
# b3lyp/6-31G Opt=(ModRedundant, loose) SCRF=(IEFPCM, solvent=water) SCF=(MaxCycle=500) NoSymmetry

B3LYP Initial Opt
                
0 1           
{xyz_output}

{ring_dihedrals}
''')
        
# Write input dihedral scan starting from a given optimized structure (mol) (DOES NOT maintain original pucker state)
def write_dihedral_scan_unconstrained(mol, state, dihedral_indices, output_path):
    xyz_output = get_xyz_from_mol(mol)
    scan_dihedral = get_scan_dihedral(dihedral_indices)
    
    with open(f'{output_path}/{state}.com', 'w') as f:
        f.write(f'''%chk={state}.chk
%nprocs=32
%mem=30GB
# b3lyp/6-31g Opt=(ModRedundant, loose) SCRF=(IEFPCM, solvent=water) SCF=(MaxCycle=500) NoSymmetry

B3LYP Dihedral Scan

0 1
{xyz_output}

{scan_dihedral}
''')
        
# Write input dihedral scan starting from a given optimized structure (mol) (DOES NOT maintain original pucker state)
def write_dihedral_scan_constrained(mol, state, dihedral_indices, output_path, ring_atoms):
    xyz_output = get_xyz_from_mol(mol)
    ring_dihedrals = get_ring_dihedral_constraints(ring_atoms)
    scan_dihedral = get_scan_dihedral(dihedral_indices)
    
    with open(f'{output_path}/{state}.com', 'w') as f:
        f.write(f'''%chk={state}.chk
%nprocs=32
%mem=30GB
# b3lyp/6-31g Opt=(ModRedundant, loose) SCRF=(IEFPCM, solvent=water) SCF=(MaxCycle=500) NoSymmetry

B3LYP Dihedral Scan

0 1
{xyz_output}

{ring_dihedrals}
{scan_dihedral}
''')

# Writes geometry optimization/freq analysis file for final dihedral-optimized geometry (maintains original pucker state)
def write_final_optimization_constrained(output_path, state, xyz_output, ring_atoms):
    lines = xyz_output.splitlines()
    output_lines = [line for i, line in enumerate(lines) if i > 1]
    output = '\n'.join(output_lines)

    if ring_atoms == 5:
        ring_dihedrals = 'D 1 2 3 4 F\nD 3 4 5 1 F'
    else:
        ring_dihedrals = 'D 1 2 3 4 F\nD 2 3 4 5 F\nD 4 5 6 1 F'
    
    with open(f'{output_path}/{state}.com', 'w') as f:
        f.write(f'''%chk={state}.chk
%nprocs=32
%mem=30GB
# MP2/6-311g Opt=(ModRedundant, tight, CalcAll) SCRF=(IEFPCM, solvent=water) SCF=(MaxCycle=500) NoSymmetry
                
Final Optimization Constrained

0 1
{output}

{ring_dihedrals}
''')

# Writes geometry optimization/freq analysis file for final dihedral-optimized geometry (DOES NOT maintain original pucker state)
def write_final_optimization_unconstrained(output_path, state, xyz_output):
    lines = xyz_output.splitlines()
    output_lines = [line for i, line in enumerate(lines) if i > 1]
    output = '\n'.join(output_lines)
    
    with open(f'{output_path}/{state}.com', 'w') as f:
        f.write(f'''%chk={state}.chk
%nprocs=32
%mem=30GB
# MP2/6-311g Opt=(ModRedundant, tight, CalcAll) SCRF=(IEFPCM, solvent=water) SCF=(MaxCycle=500) NoSymmetry
                
Final Optimization Unconstrained

0 1
{output}

''')
        
# Create the base molecule for a given pucker state
def make_base_mol(pucker_params, ring_atoms):
    # Assign total bond length and bond angles from parameters in Lucian Chan's paper
    if ring_atoms == 5:
        base_mol = Chem.MolFromSmiles("C1CCCO1")
        bond_len = [1.51]*5
        bond_angle = [1.816]*5
    else:
        base_mol = Chem.MolFromSmiles("O1CCCCC1")
        bond_len = [1.51]*6
        bond_angle = [1.911]*6

    AllChem.EmbedMolecule(base_mol, randomSeed=0xf00d)
    base_mol = AllChem.AddHs(base_mol)

    # Generate puckered base ring and store the coordinates of each atom
    # Argument form: RR.SetRingPuckerCoords(mol, ring atom idx, [amplitudes], [phase angles], bondlength, bondangle)
    if ring_atoms == 5:
        coords = RR.SetRingPuckerCoords(base_mol, list(range(ring_atoms)), [pucker_params[0]], [pucker_params[1]], bond_len, bond_angle)
    else:
        coords = RR.SetRingPuckerCoords(base_mol, list(range(ring_atoms)), [pucker_params[0], pucker_params[1]], [pucker_params[2]], bond_len, bond_angle)
    
    # Create a puckered state for the base ring
    for i in range(ring_atoms):
        base_mol.GetConformer().SetAtomPosition(i, coords[i])

    base_mol = optimize_H2(base_mol)

    phi = np.round(np.rad2deg(RA.GetRingPuckerCoords(coords)[1]))
    phi += 360 if phi<0 else 0
    Q = round(RA.GetRingPuckerCoords(coords)[0][0], 3)

    return base_mol, coords

# Optimizes structure of the initial pucker state (zero-point energy conformation)
def optimize_initial_structure(base_mol, pucker_mol, ring_atoms, state, output_path):
    mol_properties = AllChem.MMFFGetMoleculeProperties(pucker_mol)
    ff = AllChem.MMFFGetMoleculeForceField(pucker_mol, mol_properties)

    for j in range(ring_atoms):
        ff.AddFixedPoint(j)

    ff.Minimize()

    Chem.MolToPDBFile(pucker_mol, f'{output_path}/all_inputs/{state}.pdb')
    Chem.MolToPDBFile(base_mol, f'{output_path}/{state}/{state}_base_mol.pdb')
    write_gaussian_initial(pucker_mol, ring_atoms, state, output_path)

# Constructs two molecules from two unique smile strings:
# 1. Desired molecule - base ring with added functional groups (i.e. methoxyapiose)
# 2. Base five-membered ring (i.e. furanose)
def construct_pucker(pucker_params, state, mol, output_path):
    ring_atoms = int(mol.Ring_Atoms)
    base_mol, coords = make_base_mol(pucker_params, ring_atoms)
    pucker_mol, _ = renumber_sugar(mol.SMILE, ring_atoms)
    pucker_mol = Chem.AddHs(pucker_mol)
    pucker_mol = align(base_mol, pucker_mol, ring_atoms)

    atom_map = [[i,i] for i in range(ring_atoms)]
    Chem.rdMolAlign.AlignMol(pucker_mol, base_mol, atomMap=atom_map)
    
    # Assign the aligned coordinates to the desired molecule
    for i in range(ring_atoms):
        pucker_mol.GetConformer().SetAtomPosition(i, coords[i])
        
    pucker_mol = optimize_H2(pucker_mol)

    optimize_initial_structure(base_mol, pucker_mol, ring_atoms, state, output_path)

# Iterates through initial optimization folder and generates dihedral scan parameters for each rotable dihedral
def construct_first_dihedrals(state, mol, output_path, optimized_folder_path):
    ring_atoms = int(mol.Ring_Atoms)
    opt_file_path = glob.glob(f'{optimized_folder_path}/{ring_atoms}_mem/{mol.Name}/outputs/xyz_outputs/{state}*opt.xyz')
    _, atom_indices = renumber_sugar(mol.SMILE, ring_atoms)

    for i, indices in enumerate(atom_indices):
        opt_mol = get_mol_from_xyz(opt_file_path[0])
        substituent_number = i + 1
        Chem.MolToPDBFile(opt_mol, f'{output_path}/base_mol_{substituent_number}.pdb')
        new_output_path = f'{output_path}/{substituent_number}'
        os.system(f'mkdir {new_output_path}')
        write_dihedral_scan_constrained(opt_mol, state, indices, new_output_path, ring_atoms)

# Iterates through a previous dihedral scan folder and generates dihedral scan parameters for each rotable dihedral
# Excludes the dihedral that was last optimized to save on computational costs
def construct_repeat_dihedrals(mol, output_path, optimized_folder_path):
    ring_atoms = int(mol.Ring_Atoms)
    all_initial_structures = glob.glob(f'{optimized_folder_path}/{ring_atoms}_mem/{mol.Name}/outputs/Next_Scan_Starting_Points/*.xyz')
    for structure in all_initial_structures:
        name = os.path.basename(structure).split('_')
        state = name[0]
        substituent = int(name[1])
        output_state = f'{output_path}/{state}'
        os.system(f'mkdir {output_state}')
        _, atom_indices = renumber_sugar(mol.SMILE, ring_atoms)

        for i, indices in enumerate(atom_indices):
            opt_mol = get_mol_from_xyz(structure)
            substituent_number = i + 1
            if substituent != substituent_number:
                Chem.MolToPDBFile(opt_mol, f'{output_state}/base_mol_{substituent_number}.pdb')
                new_output_path = f'{output_state}/{substituent_number}'
                os.system(f'mkdir {new_output_path}')
                write_dihedral_scan_constrained(opt_mol, state, indices, new_output_path, ring_atoms)

# Copies submission scripts for the dihedral scans into the desired folder
# opt_type variable is 'initial' for an initial optimization and 'final' for a final optimization
def submission_scripts_initial_final(output_path, opt_type):
    check_5_mem = (sugar_info['Include'] == 'yes') & (sugar_info['Ring_Atoms'] == 5)
    check_6_mem = (sugar_info['Include'] == 'yes') & (sugar_info['Ring_Atoms'] == 6)

    if check_5_mem.any():
        furanose_path = os.path.join(output_path, '5_mem')

        if opt_type == 'initial':
            os.system(f'mkdir {furanose_path}')

        os.system(f'cp submission_files/{opt_type}/compile_{opt_type}_5_mem.sh {furanose_path}')
        os.system(f'cp submission_files/gaussian_submission_kestrel.sh {furanose_path}')
        os.system(f'cp submission_files/{opt_type}/run_all.sh {furanose_path}')
    
    if check_6_mem.any():
        pyranose_path = os.path.join(output_path, '6_mem')

        if opt_type == 'initial':
            os.system(f'mkdir {pyranose_path}')

        os.system(f'cp submission_files/{opt_type}/compile_{opt_type}_6_mem.sh {pyranose_path}')
        os.system(f'cp submission_files/gaussian_submission_kestrel.sh {pyranose_path}')
        os.system(f'cp submission_files/{opt_type}/run_all.sh {pyranose_path}')
        
# Copies submission scripts for the dihedral scans into the desired folder
def submission_scripts_dihedral(output_path):
    check_5_mem = (sugar_info['Include'] == 'yes') & (sugar_info['Ring_Atoms'] == 5)
    check_6_mem = (sugar_info['Include'] == 'yes') & (sugar_info['Ring_Atoms'] == 6)

    if check_5_mem.any():
        furanose_path = os.path.join(output_path, '5_mem')
        os.system(f'mkdir {furanose_path}')
        os.system(f'cp submission_files/dihedrals/compile_dihedrals.py {furanose_path}')
        os.system(f'cp submission_files/dihedrals/*.sh {furanose_path}')
        os.system(f'cp submission_files/gaussian_submission_kestrel.sh {furanose_path}')

    if check_6_mem.any():
        pyranose_path = os.path.join(output_path, '6_mem')
        os.system(f'mkdir {pyranose_path}')
        os.system(f'cp submission_files/dihedrals/compile_dihedrals.py {pyranose_path}')
        os.system(f'cp submission_files/dihedrals/*.sh {pyranose_path}')
        os.system(f'cp submission_files/gaussian_submission_kestrel.sh {pyranose_path}')
    
# Determines if the sugar if 5- or 6-membered and returns the relevant input path
def get_sugar_input_path(output_path, mol):
    if int(mol.Ring_Atoms) == 5:
        input_puckers = c5_params
        input_states = c5_states
        sugar_path = os.path.join(output_path, '5_mem', mol.Name)
    else:
        input_puckers = c6_params
        input_states = c6_states
        sugar_path = os.path.join(output_path, '6_mem', mol.Name)

    return sugar_path, input_puckers, input_states

# Generates all canonical pucker states for a molecule (constrained optimization to maintain pucker state)
def generate_puckers(output_path):
    os.system(f'mkdir {output_path}')
    submission_scripts_initial_final(output_path, 'initial')
    
    for mol in sugar_info.itertuples():
        sugar_path, input_puckers, input_states = get_sugar_input_path(output_path, mol)

        if mol.Include == 'yes':
            os.system(f'mkdir {sugar_path}')
            os.system(f'mkdir {sugar_path}/all_inputs')
            for i in range(len(input_puckers[0])):
                state = input_states[i]
                state_path = os.path.join(sugar_path, state)
                pucker = [list[i] for list in input_puckers]
                os.system(f'mkdir {state_path}')
                construct_pucker(pucker, state, mol, sugar_path)

# Generates rotamers in specified degree increments starting with a given base structure
def initial_dihedral_scan(output_path, optimized_input_path):
    os.system(f'mkdir {output_path}')
    submission_scripts_dihedral(output_path)

    for mol in sugar_info.itertuples():
        sugar_path, input_puckers, input_states = get_sugar_input_path(output_path, mol)

        if mol.Include == 'yes':
            os.system(f'mkdir {sugar_path}')
            for i in range(len(input_puckers[0])):
                state = input_states[i]
                state_path = os.path.join(sugar_path, state)
                os.system(f'mkdir {state_path}')
                construct_first_dihedrals(state, mol, state_path, optimized_input_path)

# Takes in previously optimized low-energy dihedral scan geometry and uses this as a starting point for the next dihedral scan
# Every subsequent dihedral scan following the initial scan rotates every dihedral except for the one that was previously optimized
def repeat_dihedral_scan(output_path, optimized_input_path):
    os.system(f'mkdir {output_path}')
    submission_scripts_dihedral(output_path)

    for mol in sugar_info.itertuples():
        sugar_path, _, _ = get_sugar_input_path(output_path, mol)

        if mol.Include == 'yes':
            os.system(f'mkdir {sugar_path}')
            construct_repeat_dihedrals(mol, sugar_path, optimized_input_path)

# Generates final optimizations at MP2 level
# Option 1: ring_atoms value of 5 or 6 for constrained optimization
# Option 2: ring_atoms value of 0 for unconstrained optimization
def final_opt(output_path, ring_atoms):
    submission_scripts_initial_final(output_path, 'final')

    for mol in sugar_info.itertuples():
        sugar_path, _, _ = get_sugar_input_path(output_path, mol)

        if mol.Include == 'yes':
            input_path = os.path.join(sugar_path, 'all_inputs')

            for xyz_file in os.listdir(input_path):
                state = xyz_file.split('_')[0]
                state_path = os.path.join(sugar_path, state)
                os.system(f'mkdir {state_path}')
                xyz_path = os.path.join(input_path, xyz_file)
                mol = get_mol_from_xyz(xyz_path)
                xyz_block = Chem.MolToXYZBlock(mol)

                if ring_atoms != 0:
                    write_final_optimization_constrained(state_path, state, xyz_block, ring_atoms)
                else:
                    write_final_optimization_unconstrained(state_path, state, xyz_block)

# Helper function for final optimizations at MP2 level (constrained optimization to maintain pucker state)
def final_constrained(master_path, ring_atoms):
    final_opt(master_path, ring_atoms)

# Helper function for final optimizations at MP2 level a molecule (unconstrained optimization - does not maintain pucker state)
def final_unconstrained(master_path):
    final_opt(master_path, 0)