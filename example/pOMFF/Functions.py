# -*- coding: utf-8 -*-

import re
import os
from collections import defaultdict
import subprocess
def flatten_to_strings(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):  
            flat_list.extend(flatten_to_strings(item))
        else:
            flat_list.append(str(item))  
    return flat_list

def check_and_delete_file(filename):
    
    target_file = filename
    
   
    folder_path = os.getcwd()
    file_path = os.path.join(folder_path, target_file)
    
    print(f"正在检查文件夹: {folder_path}") 
    
    if os.path.isfile(file_path):
        print(f"找到文件: {target_file}")
        print(f"文件路径: {file_path}")
        print(f"文件大小: {os.path.getsize(file_path)} 字节")
        
        while True:
            choice = input("是否要删除此文件？(y/n): ").lower()
            if choice == 'y':
                try:
                    os.remove(file_path)
                    print("文件已成功删除")
                    break
                except Exception as e:
                    print(f"删除文件时出错: {e}")
                    break
            elif choice == 'n':
                print("文件保留未删除")
                break
            else:
                print("请输入 y 或 n")
    else:
        print(f"文件夹中未找到文件: {target_file}")


# 原子质量字典
ATOMIC_MASS ={'H':1.008,'He':4.003,'Li':6.940,'Be':9.012,'B':10.810,'C':12.011,'N':14.007,'O':15.999,'F':18.998,'Ne':20.180,
'Na':22.990,'Mg':24.305,'Al':26.982,'Si':28.086,'P':30.974,'S':32.060,'Cl':35.450,'Ca':40.078,'Br':79.904,'I':126.90}

def parse_chemical_formula(formula):
    """
    eg"C2H4" -> {'C': 2, 'H': 4}
    "CH3COOH" -> {'C': 2, 'H': 4, 'O': 2}
    """
    pattern = r'([A-Z][a-z]*)(\d*)'
    elements = re.findall(pattern, formula)
    element_counts = {}
    
    for elem, count in elements:
        count = int(count) if count else 1  
        element_counts[elem] = element_counts.get(elem, 0) + count
    
    return element_counts

def calculate_molecular_mass(element_counts):
    """
    根据元素计数字典计算分子质量
    :param element_counts: 字典，格式如 {'C': 2, 'H': 4}
    :return: 分子质量（float）
    :raises: ValueError 如果包含未知元素
    
    Calculate molecular mass from an element count dictionary
    :param element_counts: Dictionary with format such as {'C': 2, 'H': 4}
    :return: Molecular mass (float)
    :raises: ValueError if an unknown element is included
    
    
    """
    total_mass = 0.0
    
    for element, count in element_counts.items():
        if element in ATOMIC_MASS:
            total_mass += ATOMIC_MASS[element] * count
        else:
            raise ValueError(f"未知元素符号: {element}")
    
    return total_mass

def calculate_mass_from_formula(formula):
    """
    直接从化学式字符串计算分子质量（整合函数）
    :param formula: 化学式字符串，如 "C2H4"
    :return: 分子质量（float）
    
    Calculate molecular mass directly from a chemical formula string (integrated function)
    :param formula: Chemical formula string, e.g., "C2H4"
    :return: Molecular mass (float)
    
    """
    element_counts = parse_chemical_formula(formula)
    #print(element_counts)
    return calculate_molecular_mass(element_counts)

def bonds2angelsdihedrals(bonds):

    adjacency = defaultdict(list)
    for a1, a2 in bonds:
        adjacency[a1].append(a2)
        adjacency[a2].append(a1)
    
    angles = []
    for a2 in adjacency:
        neighbors = adjacency[a2]
        if len(neighbors) >= 2:  
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    a1, a3 = neighbors[i], neighbors[j]
                    angles.append([a1, a2, a3])  
    
    
    dihedrals = []
    visited_dihedrals = set()
    for a2 in adjacency:
        for a3 in adjacency[a2]:
            a1_candidates = [a for a in adjacency[a2] if a != a3]
            a4_candidates = [a for a in adjacency[a3] if a != a2]
            
            if a1_candidates and a4_candidates:
                for a1 in a1_candidates:
                    for a4 in a4_candidates:
                        if a1 != a4:  
                            dihedral1 = (a1, a2, a3, a4)
                            dihedral2 = (a4, a3, a2, a1)
                            
                            if dihedral1 not in visited_dihedrals and \
                               dihedral2 not in visited_dihedrals:
                                dihedrals.append([a1, a2, a3, a4])
                                visited_dihedrals.add(dihedral1)
                                visited_dihedrals.add(dihedral2)
    
    return angles, dihedrals          
    
def load_parameters(param_files):
    """
    return: bond_params, angle_params, dihedral_params
    """
    bond_params = {}
    if os.path.exists(param_files['bonds']):
        with open(param_files['bonds'], 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    type1, type2 = sorted([parts[0], parts[1]]) 
                    bond_params[(type1, type2)] = {
                        'length': float(parts[2]),
                        'k': float(parts[3])
                    }
    
    angle_params = {}
    if os.path.exists(param_files['angles']):
        with open(param_files['angles'], 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    angle_params[(parts[0], parts[1], parts[2])] = {
                        'angle': float(parts[3]),
                        'k': float(parts[4])
                    }
    
    dihedral_params = {}
    if os.path.exists(param_files['dihedrals']):
        with open(param_files['dihedrals'], 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()

                    key = tuple(parts[:4])
                    dihedral_params[key] = {
                        'phase': float(parts[4]),
                        'k': float(parts[5]),
                        'multiplicity': int(parts[6])
                    }
    
    return bond_params, angle_params, dihedral_params

def find_parameter(atom_types, params_dict, wildcard='*'):
    """
    """
    if atom_types in params_dict:
        return params_dict[atom_types]
    
    for param_types, values in params_dict.items():
        match = True
        for i in range(len(atom_types)):
            if param_types[i] != wildcard and param_types[i] != atom_types[i]:
                match = False
                break
        if match:
            return values
    
    defaults = {
        'length': 0.15, 'k': 2500.0,  
        'angle': 120.0, 'k_angle': 100.0, 
        'phase': 180.0, 'k_dihedral': 10.0, 'multiplicity': 2  
    }
    return {k: defaults.get(k.replace('k_', ''), v) for k, v in defaults.items()}

#####################################################################################   
#                                                                                   #
#                        classes and write_functions                                #
#                                                                                   #
##################################################################################### 

class polarizable_water_ions:
    def __init__(self,  molname, atom1, atom2, atom3, mass,nrexcl=1 ):
        self.nrexcl = nrexcl
        self.molname = molname  
        self.atom1 = atom1 #{"atomtype":"WO","q":0}
        self.atom2 = atom2 #{"atomtype":"WO","q":0}
        self.atom3 = atom3 #{"atomtype":"WO","q":0}
        self.mass = mass  
        self.l0=0.110
        self.func_bond=1
        self.theta0=0.0
        self.ktheta=4.20
        self.func_angle=2
    def get_atom_type(self):
        return self.atom1.get('atomtype', None)  
    def get_dummy_atom_name1(self):
        return self.atom2.get('atomtype', None)        
    def get_dummy_atom_name2(self):          
        return self.atom3.get('atomtype', None)        
    def get_atom_name(self):
        return self.atom1.get('atomtype', None)                      
   
    
def write_itp_polarizable_water_ions(filename,mol):
    itpFile = open(filename, 'w')
    itpFile.write('[ moleculetype ] \n ')
    itpFile.write('; molname       nrexcl  \n ')
    itpFile.write('%7s %5d \n' % (mol.molname,mol.nrexcl ))
    itpFile.write('\n')
    itpFile.write('[ atoms ] \n')
    itpFile.write(';id type resnr residu atom cgnr   charge    mass       \n')
    itpFile.write('%7d %7s %7d %7s %7s %7d %7.3f %7.3f \n' % (1,mol.atom1["atomtype"],1,mol.molname,mol.atom1["atomtype"],1,mol.atom1["q"],mol.mass))
    itpFile.write('%7d %7s %7d %7s %7s %7d %7.3f %7.3f \n' % (2,mol.atom2["atomtype"],2,mol.molname,mol.atom2["atomtype"],2,mol.atom2["q"],mol.mass))
    itpFile.write('%7d %7s %7d %7s %7s %7d %7.3f %7.3f \n' % (3,mol.atom3["atomtype"],3,mol.molname,mol.atom3["atomtype"],3,mol.atom3["q"],mol.mass))
    itpFile.write('\n')
    itpFile.write('[constraints]   \n')
    itpFile.write(';  i     j   funct   length   \n')
    itpFile.write('%7d %7d %7d %7.3f   \n' % (1,2,mol.func_bond,mol.l0 ))
    itpFile.write('%7d %7d %7d %7.3f   \n' % (1,3,mol.func_bond,mol.l0 ))
    itpFile.write('\n')   
    itpFile.write('[angles]   \n')
    itpFile.write(';   i    j   k   funct  angle    fc        \n')
    itpFile.write('%7d %7d %7d %7d %7.3f %7.3f  \n' % (2,1,3,mol.func_angle,mol.theta0,mol.ktheta))
    itpFile.write('\n')
    itpFile.write('[exclusions]  \n')
    itpFile.write('%7d %7d %7d \n' % (1,2,3))
    itpFile.write('%7d %7d %7d \n' % (2,3,1))
    itpFile.write('%7d %7d %7d \n' % (3,2,1))
    itpFile.write('\n')
    itpFile.write('\n')

class normal_molecule: ################################################################################
    def __init__(self, molname, atoms, bonds=None, angles=None, dihedrals=None,nrexcl=1,posres=None):
        self.nrexcl = nrexcl
        self.molname=molname
        self.atoms=atoms #atomtype,resid,q,mass;atomname,resname, 
        self.bonds=bonds #pair,para
        self.angles=angles #pair,para
        self.dihedrals=dihedrals #pair,para
        self.func_bond=1
        self.func_angle=2
        self.func_dihedral=1
        self.posres=posres
        self.nrexcl=nrexcl

    def get_atom_type(self):
        atoms_list=[]
        for i in range(len(self.atoms)):            
            atoms_list.append(self.atoms[i].get('atomtype'))
        return atoms_list  

    def get_atom_name(self):        
        atoms_list=[]
        for i in range(len(self.atoms)):            
            atoms_list.append(self.atoms[i].get('atomname', self.atoms[i]["atomtype"]))
        return atoms_list  

def write_itp_normal_molecule(filename,mol):
    itpFile = open(filename, 'w')
    itpFile.write('[ moleculetype ] \n ')
    itpFile.write('; molname       nrexcl  \n ')
    itpFile.write('%7s %5d \n' % (mol.molname,mol.nrexcl ))
    itpFile.write('\n')
    
    itpFile.write('[ atoms ] \n')
    itpFile.write(';id type resnr residu atom cgnr   charge    mass       \n')
        
    if mol.bonds != None:
        atom_type_map={}
        for i in range(len(mol.atoms)): 
            itpFile.write('%7d %7s %7d %7s %7s %7d %7.3f %7.3f \n' % (i+1,mol.atoms[i]["atomtype"],mol.atoms[i]['resid'],
                                                                    mol.atoms[i].get('resname', mol.molname),
                                                                    mol.atoms[i].get('atomname', mol.atoms[i]["atomtype"]),
                                                                    i+1,mol.atoms[i]["q"],mol.atoms[i]['mass']
                                                                    )
                        )             
            atom_type_map[i+1] = mol.atoms[i]["atomtype"]   
    
    
        param_files = {
            'bonds': 'bond_params.dat',
            'angles': 'angle_params.dat',
            'dihedrals': 'dihedral_params.dat'
        }
        bond_params, angle_params, dihedral_params = load_parameters(param_files)
    
    
        
        bond_list=[]
        for i in range(len(mol.bonds)):
            bond_list.append(sorted([mol.bonds[i][0],mol.bonds[i][1]]))
        #print('bond_list:',bond_list)
        #print('bonds:',bonds)          
        
        itpFile.write('[ bonds ] \n')
        itpFile.write('; i  j  funct  rb0    kb\n')     
        for a1, a2 in bond_list:
            type1, type2 = sorted([atom_type_map[a1], atom_type_map[a2]])
            list1=sorted([a1,a2]) 
            params={}        
            try:
                params['length']=mol.bonds[bond_list.index(list1)][2]
                params['k']=mol.bonds[bond_list.index(list1)][3]        
            except Exception:   
                params = find_parameter((type1, type2), bond_params)
            itpFile.write('%7d %7d %7d %7.3f %7.0f \n'%(a1,a2,mol.func_bond,params['length'],params['k']))

        angle_list_all=bonds2angelsdihedrals(bond_list)[0]
        
        itpFile.write('[ angles ] \n')
        itpFile.write(';  i    j    k   funct  angle     fc. \n')
        
        if mol.angles != None:
            angle_list=[]
            for i in range(len(mol.angles)):
                try:
                    angle_list.append(sorted([mol.angles[i][0],mol.angles[i][1],mol.angles[i][2]]))
                except Exception:
                    aa=1
            print('angle_list_all:',angle_list_all)
            
            
            for i in range(len(angle_list_all)):
                a1=angle_list_all[i][0]
                a2=angle_list_all[i][1]
                a3=angle_list_all[i][2]
                type1 = atom_type_map[a1]
                type2 = atom_type_map[a2]
                type3 = atom_type_map[a3]
                #print(type1,type2,type3)
                list1=sorted([a1,a2,a3]) 
                params={}
                if list1 in angle_list:
                    try:
                        params['angle']=mol.angles[angle_list.index(list1)][3]
                        params['k']=mol.angles[angle_list.index(list1)][4]
                        itpFile.write('%7d %7d %7d %7d %7.3f %7.0f \n'%(a1,a2,a3,mol.func_angle,params['angle'],params['k']))
                        print('special constraints for:',list1,'in',mol.molname)
                    except Exception:
                        print('no constraints for:',list1,'in',mol.molname)                      
                else:
                    params = find_parameter((type1, type2, type3), angle_params)
                    itpFile.write('%7d %7d %7d %7d %7.3f %7.0f \n'%(a1,a2,a3,mol.func_angle,params['angle'],params['k']))           
            
        else:
            for i in range(len(angle_list_all)):
                a1=angle_list_all[i][0]
                a2=angle_list_all[i][1]
                a3=angle_list_all[i][2]
                type1 = atom_type_map[a1]
                type2 = atom_type_map[a2]
                type3 = atom_type_map[a3]
                params = find_parameter((type1, type2, type3), angle_params)
                itpFile.write('%7d %7d %7d %7d %7.3f %7.0f \n'%(a1,a2,a3,mol.func_angle,params['angle'],params['k']))
            
        dihedral_list_all=bonds2angelsdihedrals(bond_list)[1]
        print(dihedral_list_all)
        itpFile.write('[ dihedrals ] \n')
        itpFile.write('; ai	aj	ak	al	funct	phi0	cp	mult \n')
        
        if mol.dihedrals != None :
            if mol.dihedrals == 'NoConstraints':
                print('No dihedral constraint')
            else:
                dihedral_list=[]
                for i in range(len(mol.dihedrals)):
                    try:
                        dihedral_list.append(sorted([mol.dihedrals[i][0],mol.dihedrals[i][1],mol.dihedrals[i][2],mol.dihedrals[i][3]]))
                    except Exception:
                        aa=1
                #print('dihedral_list:',dihedral_list)
                #print('dihedrals:',dihedrals)
            
                for i in range(len(dihedral_list_all)):
                    a1=dihedral_list_all[i][0]
                    a2=dihedral_list_all[i][1]
                    a3=dihedral_list_all[i][2]
                    a4=dihedral_list_all[i][3]
                    type1 = atom_type_map[a1]
                    type2 = atom_type_map[a2]
                    type3 = atom_type_map[a3]
                    type4 = atom_type_map[a4]
                                
                    list1=sorted([a1,a2,a3,a4])
                    params={}
                    if list1 in dihedral_list:
                        try:
                            params['phase']=mol.dihedrals[dihedral_list.index(list1)][4]
                            params['k']=mol.dihedrals[dihedral_list.index(list1)][5]
                            params['multiplicity']=mol.dihedrals[dihedral_list.index(list1)][6]
                            itpFile.write('%7d %7d %7d %7d %7d %7.3f %7.0f %7d \n'%(a1,a2,a3,a4,mol.func_dihedral,params['phase'],params['k'],params['multiplicity']))
                            print('special constraints for:',list1,'in',mol.molname)
                        except Exception:
                            print('no constraints for:',list1,'in',mol.molname)                         
                    else:
                        params = find_parameter((type1, type2, type3, type4),  dihedral_params)
                        itpFile.write('%7d %7d %7d %7d %7d %7.3f %7.0f %7d \n'%(a1,a2,a3,a4,mol.func_dihedral,params['phase'],params['k'],params['multiplicity']))
        else:
            for i in range(len(dihedral_list_all)):
                a1=dihedral_list_all[i][0]
                a2=dihedral_list_all[i][1]
                a3=dihedral_list_all[i][2]
                a4=dihedral_list_all[i][3]
                type1 = atom_type_map[a1]
                type2 = atom_type_map[a2]
                type3 = atom_type_map[a3]
                type4 = atom_type_map[a4]
                params = find_parameter((type1, type2, type3, type4),  dihedral_params)
                itpFile.write('%7d %7d %7d %7d %7d %7.3f %7.0f %7d \n'%(a1,a2,a3,a4,mol.func_dihedral,params['phase'],params['k'],params['multiplicity']))
            

    else:
        itpFile.write('%7d %7s %7d %7s %7s %7d %7.3f %7.3f \n' % (1,mol.atoms[0]["atomtype"],mol.atoms[0]['resid'],
                                                                    mol.atoms[0].get('resname', mol.molname),
                                                                    mol.atoms[0].get('atomname', mol.atoms[0]["atomtype"]),
                                                                    1,mol.atoms[0]["q"],mol.atoms[0]['mass']
                                                                  )
                      )
    if mol.posres != None:
        itpFile.write('\n \n #ifdef POSRES \n')
        itpFile.write('[ position_restraints ] \n')
        itpFile.write(';ai funct fc \n')
        for i in range(len(mol.posres)):
            itpFile.write('%7d %40s \n' % (mol.posres[i],'1 0.0 0.0 POSRES_FC_LIPID'))
        itpFile.write('#endif \n')
    
    
def write_topol_top(molecules,systemname,filename="topol.top"):
    topFile = open(filename,'w')
    topFile.write(';;Generated by pOMFF.py which is writen by Yi Dong \n') #################################################
    topFile.write(';;pOMFF force field \n')
    topFile.write(';;Correspondance:lhgao@bnu.edu.cn \n')
    topFile.write('#include "pOMFF.ff/pOMFF.itp" \n')
    for i in range(len(molecules)):
        topFile.write('#include "pOMFF.ff/%s.itp"\n' % molecules[i][0].molname)  ##########################
    
    topFile.write('\n[ system ]\n')
    topFile.write('; Name\n')
    topFile.write('%10s\n' % (systemname))
    
    topFile.write('\n[ molecules ]\n')
    topFile.write('; Compound\t#mols\n')
    for i in range(len(molecules)):
        topFile.write('%-6s\t%12d\n' % (molecules[i][0].molname,molecules[i][1]))        ########################## 


def write_ffnonbonded_itp(all_atom_types,all_dummy_atom_types):
    func_nb=1
    func_nb_para_1=0.5
    func_nb_para_2=0.5
    
    itpFile = open("ffnonbonded.itp", 'w')
    
    itpFile.write('[ atomtypes ] \n')
    itpFile.write(';name mass   charge  ptype C   A  \n')
    for i in range(len(all_atom_types)):
        itpFile.write('%7s %7.3f %7.3f %7s %7.3f %7.3f \n' % (all_atom_types[i],0.0,0.0,"A",0.0,0.0))
    itpFile.write('\n')
    for i in range(len(all_dummy_atom_types)):
        itpFile.write('%7s %7.3f %7.3f %7s %7.3f %7.3f \n' % (all_dummy_atom_types[i],0.0,0.0,"A",0.0,0.0))
    itpFile.write('\n')
    
    itpFile.write('[ nonbond_params ] \n')
    itpFile.write(';i   j  func C A \n')
    for i in range(len(all_atom_types)):
        for j in range(i,len(all_atom_types)):
            #print(i,j)
            itpFile.write('%7s %7s %7d %7.1f %7.1f \n'%(all_atom_types[i],all_atom_types[j],func_nb,func_nb_para_1,func_nb_para_2))
    itpFile.write('\n')
    

def write_main_itp():
    itpFile = open("pOMFF.itp", 'w')
    
    itpFile.write('[ defaults ] \n')
    itpFile.write(';nbfunc   comb-rule   gen-pairs   fudgeLJ   fudgeQQ \n')
    itpFile.write('%7d %7d %7s %7.1f %7.1f \n' % (1, 1, "no", 1.0, 1.0))
    itpFile.write('\n')
    itpFile.write('\n')
    itpFile.write('#include "./ffnonbonded.itp" \n ')

def write_em_mdp(filename,all_atom_types,all_dummy_atom_types,emtol,tau_t,Temperature):
    mdpFile = open(filename,'w')
    mdpFile.write('define = \n')##
    mdpFile.write('title                    = pOMFF \n')
    mdpFile.write('integrator               = steep\n')
    mdpFile.write('emtol                    = %d\n' % (emtol)) #%
    mdpFile.write('dt                       = 0.0001\n')
    mdpFile.write('nsteps                   = 200000\n')
    mdpFile.write('nstcomm                  = 100\n')
    mdpFile.write('nstxout                  = 0\n')
    mdpFile.write('nstvout                  = 0\n')
    mdpFile.write('nstfout                  = 0\n')
    mdpFile.write('nstlog                   = 1000\n')
    mdpFile.write('nstenergy                = 1000\n')
    mdpFile.write('nstxtcout                = 1000\n')
    mdpFile.write('xtc_precision            = 1000\n')
    all_atom_types_dummy_atoms=all_atom_types+all_dummy_atom_types
    energygrps=''
    for i in all_atom_types_dummy_atoms:
        energygrps=energygrps+' '+i+ ' '
    mdpFile.write('energygrps               = %s\n' % (energygrps))
    energygrp_table=''
    for i in range(len(all_atom_types)):
        for j in range(i,len(all_atom_types)):
            energygrp_table=energygrp_table+all_atom_types[i]+' '+all_atom_types[j]+' '
    mdpFile.write('energygrp_table          = %s\n' % (energygrp_table))
    mdpFile.write('cutoff-scheme            = group\n') ##
    mdpFile.write('nstlist                  = 10\n')
    mdpFile.write('ns_type                  = grid\n')
    mdpFile.write('pbc                      = xyz\n')
    mdpFile.write('rlist                    = 1.35\n')
    mdpFile.write('coulombtype              = PME\n')
    mdpFile.write('rcoulomb                 = 1.35\n')
    mdpFile.write('epsilon_r                = 3.2\n')
    mdpFile.write('epsilon_rf               = 0\n')
    mdpFile.write('vdw_type                 = user\n')
    mdpFile.write('vdw-modifier             = Potential-shift-verlet\n') ##
    mdpFile.write('rvdw                     = 1.35\n')
    mdpFile.write('tcoupl                   = V-rescale \n')
    mdpFile.write('tc-grps                  = %s \n' %(energygrps)) 
    tau_t_str=''
    ref_t=''
    for i in range(len(all_atom_types_dummy_atoms)):
        tau_t_str=tau_t_str+' '+str(tau_t)+' ' 
        ref_t=ref_t+' '+str(Temperature)+' ' 
    mdpFile.write('tau_t                    = %s \n' % (tau_t_str))
    mdpFile.write('ref_t                    = %s \n' % (ref_t))
    mdpFile.write('Pcoupl                   = Parrinello-Rahman \n')
    mdpFile.write('Pcoupltype               = isotropic\n')
    mdpFile.write('tau_p                    = 12.0\n')
    mdpFile.write('compressibility          = 3e-4\n')
    mdpFile.write('ref_p                    = 1.0\n')
    mdpFile.write('gen_vel                  = no\n')
    mdpFile.write('gen_temp                 = 320\n')
    mdpFile.write('gen_seed                 = 183511\n')
    mdpFile.write('constraints              = none\n')
    mdpFile.write('constraint_algorithm     = Lincs\n')


def write_md_mdp(filename,all_atom_types,all_dummy_atom_types,dt,simu_time,tau_t,Temperature,pcouple,pcoupletype,define=None,tc_grps=None,tension=None): 
    mdpFile = open(filename,'w')
    if define == None : 
        mdpFile.write('define = \n')
    else :
        mdpFile.write('define = %s \n' % (define))
    mdpFile.write('title                    = pOMFF \n')
    mdpFile.write('integrator               = sd\n')
    mdpFile.write('dt                       = %f\n' % (dt))
    nsteps=int(simu_time/dt)
    mdpFile.write('nsteps                   = %d\n' % (nsteps)) ########
    mdpFile.write('nstcomm                  = 1000\n')
    mdpFile.write('nstxout                  = 0\n')
    mdpFile.write('nstvout                  = 0\n')
    mdpFile.write('nstfout                  = 0\n')
    mdpFile.write('nstlog                   = 1000\n')
    mdpFile.write('nstenergy                = 1000\n')
    mdpFile.write('nstxtcout                = 1000\n')
    mdpFile.write('xtc_precision            = 1000\n')
    all_atom_types_dummy_atoms=all_atom_types+all_dummy_atom_types
    energygrps=''
    for i in all_atom_types_dummy_atoms:
        energygrps=energygrps+' '+i+ ' '
    mdpFile.write('energygrps               = %s\n' % (energygrps))
    energygrp_table=''
    for i in range(len(all_atom_types)):
        for j in range(i,len(all_atom_types)):
            energygrp_table=energygrp_table+all_atom_types[i]+' '+all_atom_types[j]+' '
    mdpFile.write('energygrp_table          = %s\n' % (energygrp_table))
    mdpFile.write('cutoff-scheme            = group\n') ##
    mdpFile.write('nstlist                  = 10\n')
    mdpFile.write('ns_type                  = grid\n')
    mdpFile.write('pbc                      = xyz\n')
    mdpFile.write('rlist                    = 1.35\n')
    mdpFile.write('coulombtype              = PME\n')
    mdpFile.write('rcoulomb                 = 1.35\n')
    mdpFile.write('epsilon_r                = 3.2\n')
    mdpFile.write('epsilon_rf               = 0\n')
    mdpFile.write('vdw_type                 = user\n')
    mdpFile.write('vdw-modifier             = Potential-shift-verlet\n') ##
    mdpFile.write('rvdw                     = 1.35\n')
    mdpFile.write('tcoupl                   = V-rescale \n')
    if tc_grps == None :
        mdpFile.write('tc-grps                  = %s \n' %(energygrps)) 
        tau_t_str=''
        ref_t=''
        for i in range(len(all_atom_types_dummy_atoms)):
            tau_t_str=tau_t_str+' '+str(tau_t)+' ' 
            ref_t=ref_t+' '+str(Temperature)+' ' 
        mdpFile.write('tau_t                    = %s \n' % (tau_t_str))
        mdpFile.write('ref_t                    = %s \n' % (ref_t))
    else:
        tc_grps_ss=''
        for i in range(len(tc_grps)):
            tc_grps_ss=tc_grps_ss+tc_grps[i]+' '     
        mdpFile.write('tc-grps                  = %s \n' %(tc_grps_ss)) 
        tau_t_str=''
        ref_t=''
        for i in range(len(tc_grps)):
            tau_t_str=tau_t_str+' '+str(tau_t)+' ' 
            ref_t=ref_t+' '+str(Temperature)+' ' 
        mdpFile.write('tau_t                    = %s \n' % (tau_t_str))
        mdpFile.write('ref_t                    = %s \n' % (ref_t))        
    
    if pcouple.lower()!="no".lower():
        mdpFile.write('Pcoupl                   = %s \n' % (pcouple))
        mdpFile.write('Pcoupltype               = %s\n' % (pcoupletype))

        if pcouple.lower() == 'Parrinello-Rahman'.lower()  :    
            tau_p=12.0
        elif pcouple.lower() == 'berendsen'.lower() :
            tau_p=5.0        
        mdpFile.write('tau_p                    = %d\n' % (tau_p))
        
        if pcoupletype.lower()=='isotropic'.lower() :
            compressibility='3e-4'
            ref_p='1.0'
        elif pcoupletype.lower()=='semiisotropic'.lower():
            compressibility='3e-4  3e-4'
            ref_p='1.0  1.0 '  
        elif pcoupletype.lower()=='surface-tension'.lower():
            compressibility='3e-4 0'
            tension=tension*20
            ref_p=str(tension)+'  0.0 '  
            
        mdpFile.write('compressibility          = %s\n' % (compressibility))
        mdpFile.write('ref_p                    = %s\n' % (ref_p))
    else:
        mdpFile.write('Pcoupl                   = %s \n' % (pcouple))
    mdpFile.write('gen_vel                  = no\n')
    mdpFile.write('gen_temp                 = 320\n')
    mdpFile.write('gen_seed                 = 183511\n')
    mdpFile.write('constraints              = none\n')
    mdpFile.write('constraint_algorithm     = Lincs\n')
    if dt == 0.02:
        mdpFile.write('lincs_order              = 4\n')
        mdpFile.write('lincs_warnangle          = 30\n')



def write_make_tablessh(all_atom_types):
    filename="make_tables.sh"
    str=''
    for i in range(len(all_atom_types)):
        str=str+all_atom_types[i]+' '
        

    maketableFile= open(filename,'w')
    maketableFile.write('#!/bin/bash \n')

    maketableFile.write('chmod +x ./TableA; \n')
    maketableFile.write('chmod +x ./TableM; \n \n')
    maketableFile.write('''	echo 'nam1 nam2 ep R0 al be lambda '>table_list_use.dat \n''')
    maketableFile.write('list=(%s) \n \n' % (str))

    maketableFile.write('''n_tt=$(awk 'END{print NR}' table_list.dat) \n \n''')

    maketableFile.write('''for tt in $(seq 3 1 $n_tt) \n''')
    maketableFile.write('''do \n''')
    maketableFile.write('''table=$(cat ./table_list.dat | awk NR==$tt'{ print $1 }') \n''')
    maketableFile.write('''nam1=$(cat ./table_list.dat | awk NR==$tt'{ print $3 }') \n''')
    maketableFile.write('''nam2=$(cat ./table_list.dat | awk NR==$tt'{ print $4 }') \n''')
    maketableFile.write('''tep=$(cat ./table_list.dat | awk NR==$tt'{ print $5 }') \n''')
    maketableFile.write('''tR0=$(cat ./table_list.dat | awk NR==$tt'{ print $6 }') \n''')
    maketableFile.write('''tal=$(cat ./table_list.dat | awk NR==$tt'{ print $7 }') \n''')
    maketableFile.write('''tbe=$(cat ./table_list.dat | awk NR==$tt'{ print $8 }') \n''')
    maketableFile.write('''lambda=$(cat ./table_list.dat | awk NR==$tt'{ print $9 }') \n \n''')

    maketableFile.write('''if [[ "${list[@]}"  =~ "${nam1}" ]] && [[ "${list[@]}"  =~ "${nam2}" ]] ; then \n''')
    maketableFile.write('''	lists=($nam1 $nam2) \n''')
    maketableFile.write('''	sorted_lists=($(printf "%s\\n" "${lists[@]}" | sort -f)) \n''')
    maketableFile.write('''	echo $nam1 $nam2 $tep $tR0 $tal $tbe $lambda >>table_list_use.dat \n''')
    maketableFile.write('''	if [ $table == 1 ] ;then \n''')
    maketableFile.write('''		echo -e "$tep $tR0 $tal $tbe $lambda 1.35" |./TableA; mv table.xvg table_${sorted_lists[0]}_${sorted_lists[1]}.xvg \n''')
#    maketableFile.write('''		mv tabled.xvg tabled_${sorted_lists[0]}_${sorted_lists[1]}.xvg \n''')
    maketableFile.write('''		echo $tt $nam1 $nam2 \n''')
    maketableFile.write('''	elif [ $table == 2 ];then \n''')
    maketableFile.write('''		echo -e "$tep $tR0 $tal $tbe $lambda  1.35" |./TableM; mv table.xvg table_${sorted_lists[0]}_${sorted_lists[1]}.xvg \n''')
#    maketableFile.write('''		mv tabled.xvg tabled_${sorted_lists[0]}_${sorted_lists[1]}.xvg \n''')
    maketableFile.write('''		echo $tt  $nam1 $nam2 \n''')
    maketableFile.write('''	fi \n''')
    maketableFile.write('''fi \n''')
    maketableFile.write('''done \n \n''')
    maketableFile.write('''echo -e "0.000 0.000 0.000 0.000 1.0 1.35" |./TableM; \n''')
    



    

def find_atom_types_names(molecules,polarizable_molecules):
    if polarizable_molecules != None: 
    
        all_atom_types=[]
        for i in range(len(molecules)):
            all_atom_types.append(molecules[i][0].get_atom_type())
        all_dummy_atom_types = [s.get_dummy_atom_name1() for s in polarizable_molecules if s.get_dummy_atom_name1() is not None] +\
                            [s.get_dummy_atom_name2() for s in polarizable_molecules if s.get_dummy_atom_name2() is not None]
        
        all_atom_types = list(set(flatten_to_strings(all_atom_types)))
        all_atom_types.sort(key=lambda x: x.lower())
        all_dummy_atom_types = list(set(all_dummy_atom_types))
        all_dummy_atom_types.sort(key=lambda x: x.lower())
    
    
        #找到所有珠子名字并去重 Find all bead names and remove duplicates
        all_atom_names=[]
        for i in range(len(molecules)):
            all_atom_names.append(molecules[i][0].get_atom_name())
        all_atom_names = list(set(flatten_to_strings(all_atom_names)))
        all_atom_names.extend(all_dummy_atom_types)
        all_atom_names.sort(key=lambda x: x.lower())
    else :
        all_atom_types=[]
        all_dummy_atom_types=[]
        for i in range(len(molecules)):
            all_atom_types.append(molecules[i][0].get_atom_type())
        
        all_atom_types = list(set(flatten_to_strings(all_atom_types)))
        all_atom_types.sort(key=lambda x: x.lower())
    
    
        #找到所有珠子名字并去重 Find all bead names and remove duplicates
        all_atom_names=[]
        for i in range(len(molecules)):
            all_atom_names.append(molecules[i][0].get_atom_name())
        all_atom_names = list(set(flatten_to_strings(all_atom_names)))
        all_atom_names.sort(key=lambda x: x.lower())
        
    return all_atom_types,all_dummy_atom_types,all_atom_names

def gen_files(molecules,polarizable_molecules,systemname):
    print('gen top and itp files')
    write_topol_top(molecules,systemname)
    
    #找到所有珠子类型并去重 Find all bead types and remove duplicates
        
    all_atom_types=find_atom_types_names(molecules,polarizable_molecules)[0]
    all_dummy_atom_types=find_atom_types_names(molecules,polarizable_molecules)[1]
    all_atom_names=find_atom_types_names(molecules,polarizable_molecules)[2]
 
    
    indexinpFile=open ('index.inp','w')
    for i in all_atom_names:
        indexinpFile.write('a %7s \n' % (i))
    indexinpFile.write('\n q \n \n')     
    print('################################################################')
    print("all_atom_types:", all_atom_types)  
    print("all_dummy_atom_types:", all_dummy_atom_types)    
    print('all_atom_names:',all_atom_names)
    print('################################################################')

    # 创建力场文件目录并进入
    fffolder_path = "pOMFF"
    
    try:
        os.mkdir(fffolder_path)
        print(f"文件夹 '{fffolder_path}' 创建成功")
    except FileExistsError:
        print(f"文件夹 '{fffolder_path}' 已存在")
    except Exception as e:
        print(f"创建失败: {e}")
    
    if os.path.isdir(fffolder_path):
        os.chdir(fffolder_path)  # 进入文件夹
        print(f"当前工作目录已更改为: {os.getcwd()}")
    else:
        print(f"文件夹不存在: {folder_path}")
    
    #写主文件和非键itp;forcefield.itp and ffnonbonded.itp 
    write_main_itp()
    write_ffnonbonded_itp(all_atom_types,all_dummy_atom_types)
    
    #写分子itp;molecular.itp 
    #check_and_delete_file("ions.itp")
    if polarizable_molecules != None :
        for i in range(len(polarizable_molecules)):
            filename=polarizable_molecules[i].molname+".itp"
            write_itp_polarizable_water_ions(filename,polarizable_molecules[i])
        
        for i in range(len(molecules)):
            if molecules[i][0] not in  polarizable_molecules:
                filename=molecules[i][0].molname+".itp"
                write_itp_normal_molecule(filename,molecules[i][0])
    else:
        for i in range(len(molecules)):
            filename=molecules[i][0].molname+".itp"
            write_itp_normal_molecule(filename,molecules[i][0])        
    os.chdir("..")
    print(f"当前工作目录已更改为: {os.getcwd()}") 
    
    os.system("mkdir pOMFF.ff")
    
    #fffolder_path=pOMFF
    os.system("mv pOMFF/*.itp pOMFF.ff")
    
    
    #写make_tables.sh 
    write_make_tablessh(all_atom_types)
    os.system("chmod +x ./make_tables.sh ")
