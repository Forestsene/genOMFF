# -*- coding: utf-8 -*-

import re
import os
from collections import defaultdict
import subprocess

from pOMFF.Functions import gen_files,find_atom_types_names,write_md_mdp,write_em_mdp
                       
#####################################################################################   
#                                                                                   #
#                         gen files for specific systems                            #
#                                                                                   #
#####################################################################################   

#体系中包含的分子和个数
from pOMFF.moleculars import NA,CL,PW
moleculars=([NA,9],[CL,9],[PW,400])#得写两个Write at least two.
polarizable_moleculars=[NA,CL,PW]

systemname="NaCl"

gen_files(moleculars,polarizable_moleculars,systemname)

#print('gen_mdp_files')
#找到所有珠子类型并去重

all_atom_types=find_atom_types_names(moleculars,polarizable_moleculars)[0]
all_dummy_atom_types=find_atom_types_names(moleculars,polarizable_moleculars)[1]
#print('################################################################')
#print("all_atom_types:", all_atom_types)  
#print("all_dummy_atom_types:", all_dummy_atom_types) 
#print('################################################################')

#filename,all_atom_types,all_dummy_atom_types,emtol,tau_t,Temperature
write_em_mdp('em.mdp',all_atom_types,all_dummy_atom_types,100,0.3,298)   
#filename,all_atom_types,all_dummy_atom_types,dt,simu_time,tau_t,Temperature,pcouple,pcoupletype,define=None,tc_grps=None,tension=None
write_md_mdp('npt.mdp',all_atom_types,all_dummy_atom_types,0.02,4000,0.3,298,'berendsen','isotropic')
write_md_mdp('100ns_npt.mdp',all_atom_types,all_dummy_atom_types,0.02,100000,0.3,298,'berendsen','isotropic')
write_md_mdp('100ns_npt1.mdp',all_atom_types,all_dummy_atom_types,0.02,100000,0.3,298,'berendsen','isotropic')
write_md_mdp('100ns_npt2.mdp',all_atom_types,all_dummy_atom_types,0.02,100000,0.3,298,'Parrinello-Rahman','isotropic')
write_md_mdp('nvt.mdp',all_atom_types,all_dummy_atom_types,0.02,40000,0.3,298,'no','isotropic')
write_md_mdp('100ns_nvt.mdp',all_atom_types,all_dummy_atom_types,0.02,100000,0.3,298,'no','isotropic')

write_md_mdp('1_nvt.mdp',all_atom_types,all_dummy_atom_types,0.01,1000,0.3,298,'no','isotropic',define='-DPOSRES -DPOSRES_FC_LIPID=1000',tc_grps=['SOLV','MOL'])
write_md_mdp('2_nvt.mdp',all_atom_types,all_dummy_atom_types,0.01,1000,0.3,298,'no','isotropic',define='-DPOSRES -DPOSRES_FC_LIPID=800',tc_grps=['SOLV','MOL'])
write_md_mdp('3_nvt.mdp',all_atom_types,all_dummy_atom_types,0.01,2000,0.3,298,'no','isotropic',define='-DPOSRES -DPOSRES_FC_LIPID=500',tc_grps=['SOLV','MOL'])
write_md_mdp('4_nvt.mdp',all_atom_types,all_dummy_atom_types,0.01,2000,0.3,298,'no','isotropic',define='-DPOSRES -DPOSRES_FC_LIPID=200',tc_grps=['SOLV','MOL'])
write_md_mdp('5_nvt.mdp',all_atom_types,all_dummy_atom_types,0.01,10000,0.3,298,'no','isotropic',tc_grps=['SOLV','MOL'])
write_md_mdp('md.mdp',all_atom_types,all_dummy_atom_types,0.01,200000,0.3,298,'no','isotropic',tc_grps=['SOLV','MOL'])



