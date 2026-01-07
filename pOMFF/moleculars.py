# -*- coding: utf-8 -*-

import re
import os
from collections import defaultdict
import subprocess
from .Functions import polarizable_water_ions
from .Functions import normal_molecule 
from .Functions import calculate_mass_from_formula
#####################################################################################   
#                                                                                   #
#                         molecule informations                                     #
#                                                                                   #
#####################################################################################   
#molname, atom1, atom2, atom3, mass,nrexcl=1                
PW = polarizable_water_ions("PW", {"atomtype":"WO","q":0}, {"atomtype":"WP","q":0.621}, {"atomtype":"WN","q":-0.621}, calculate_mass_from_formula("H8O4")/3 )
NA = polarizable_water_ions("NA", {"atomtype":"NAC","q":1}, {"atomtype":"NAP","q":0.5}, {"atomtype":"NAN","q":-0.5}, calculate_mass_from_formula("NaH6O3")/3 )
CL = polarizable_water_ions("CL", {"atomtype":"CLC","q":-1}, {"atomtype":"CLP","q":0.8}, {"atomtype":"CLN","q":-0.8}, calculate_mass_from_formula("ClH6O3")/3 )
BR = polarizable_water_ions("BR", {"atomtype":"BRC","q":-1}, {"atomtype":"BRP","q":0.8}, {"atomtype":"BRN","q":-0.8}, calculate_mass_from_formula("BrH6O3")/3 )


#molname, atoms, bonds=None, angles=None, dihedrals=None,nrexcl=1
#未定义的resname=molname;未定义的atomname=atomtype 
#所用成键的对都要列出，数字从小到大，列出参数的用这里的参数,没有列出的用默认参数
#angel数字从小到大，列出的参数用这里的参数，只列了对但没列参数的就没有约束，没有列出的用默认参数，None全部用默认参数
#dihedral数字从小到大，列出的参数用这里的参数，只列了对但没列参数的就没有约束,没有列出的用默认参数，None全部用默认参数,NoConstraints全部没有约束

#molname, atoms, bonds=None, angles=None, dihedrals=None, nrexcl=1
#Undefined residue name = molname; undefined atom name = atom type.
#All bonding pairs must be listed, in ascending numerical order. Listed parameters use the specified values; unlisted ones use default parameters.
#Angles in ascending numerical order. Listed parameters use the specified values; pairs listed without parameters are unconstrained; unlisted angles use default parameters; None means all use default parameters.
#Dihedrals in ascending numerical order. Listed parameters use the specified values; pairs listed without parameters are unconstrained; unlisted dihedrals use default parameters; None means all use default parameters; NoConstraints means all are unconstrained.

MSO4 = normal_molecule("MSO4",[{"atomtype":"MSO4",'resid':1,"q":-1,'mass':calculate_mass_from_formula("CH3SO4")}])
MSO4H =normal_molecule("MSO4H",[{"atomtype":"MSO4H",'resid':1,"q":0,'mass':calculate_mass_from_formula("CH3SO4H")}])
NC4 = normal_molecule("NC4",[{"atomtype":"NC4",'resid':1,"q":1,'mass':calculate_mass_from_formula("C4H12N")}])
C4 = normal_molecule("C4",[{"atomtype":"C4",'resid':1,"q":0,'mass':calculate_mass_from_formula("C4H10")}])
C3 = normal_molecule("C3",[{"atomtype":"C3",'resid':1,"q":0,'mass':calculate_mass_from_formula("C3H8")}])


SDS=normal_molecule(molname='SDS',atoms=[{'atomtype':'MSO4','resid':1,'q':-1,'mass':calculate_mass_from_formula('CH2SO4')},
                                 {'atomtype':'C04','resid':1,'q':0,'mass':calculate_mass_from_formula('C4H8')},
                                 {'atomtype':'C04','resid':1,'q':0,'mass':calculate_mass_from_formula('C4H8')},
                                 {'atomtype':'C03','resid':1,'q':0,'mass':calculate_mass_from_formula('C3H7')}
                                ],
                          bonds=[[1,2],[2,3],[3,4,0.430,2000]],
                          angles=None,
                          dihedrals=[[1,2,3,4]],
                          posres=[1]
                    )

SDS344=normal_molecule(molname='SDS344',atoms=[{'atomtype':'MSO4','resid':1,'resname':'SDS','q':-1,'mass':calculate_mass_from_formula('CH2SO4')},
                                 {'atomtype':'C03','resid':1,'resname':'SDS','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                 {'atomtype':'C04','resid':1,'resname':'SDS','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                 {'atomtype':'C04','resid':1,'resname':'SDS','q':0,'mass':calculate_mass_from_formula('C4H9')},
                                ],
                                bonds=[[1,2],[2,3,0.430,2000],[3,4]],
                                angles=None,
                                dihedrals=[[1,2,3,4]],
                                posres=[1]
                        )

MSO4C12=normal_molecule(molname='MSO4C12',atoms=[{'atomtype':'MSO4','resid':1,'resname':'MSO4C12','q':-1,'mass':calculate_mass_from_formula('CH2SO4')},
                                 {'atomtype':'C03','resid':1,'resname':'MSO4C12','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                 {'atomtype':'C03','resid':1,'resname':'MSO4C12','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                 {'atomtype':'C03','resid':1,'resname':'MSO4C12','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                 {'atomtype':'C03','resid':1,'resname':'MSO4C12','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                ],
                                bonds=[[1,2],[2,3],[3,4],[4,5]],
                                angles=None,
                                dihedrals=[[1,2,3,4]],
                                posres=[1]
                        )

DTA=normal_molecule(molname='DTA',atoms=[{'atomtype':'NC4','resid':1,'q':1,'mass':calculate_mass_from_formula('C4H11N')},
                                 {'atomtype':'C04','resid':1,'q':0,'mass':calculate_mass_from_formula('C4H8')},
                                 {'atomtype':'C04','resid':1,'q':0,'mass':calculate_mass_from_formula('C4H8')},
                                 {'atomtype':'C03','resid':1,'q':0,'mass':calculate_mass_from_formula('C3H7')}
                                ],
                           bonds=[[1,2],[2,3],[3,4]],
                           angles=None,
                           dihedrals=[[1,2,3,4]],
                           posres=[1]
                    )


DTA344=normal_molecule(molname='DTA344',atoms=[{'atomtype':'NC4','resid':1,'resname':'DTA','q':1,'mass':calculate_mass_from_formula('C4H11N')},
                                       {'atomtype':'C03','resid':1,'resname':'DTA','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C04','resid':1,'resname':'DTA','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'DTA','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3,0.430,2000],[3,4]],
                                angles=None,
                                dihedrals=[[1,2,3,4]],
                                posres=[1]
                      )

NC4C12=normal_molecule(molname='NC4C12',atoms=[{'atomtype':'NC4','resid':1,'resname':'NC4C12','q':1,'mass':calculate_mass_from_formula('C4H11N')},
                                 {'atomtype':'C03','resid':1,'resname':'NC4C12','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                 {'atomtype':'C03','resid':1,'resname':'NC4C12','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                 {'atomtype':'C03','resid':1,'resname':'NC4C12','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                 {'atomtype':'C03','resid':1,'resname':'NC4C12','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                ],
                                bonds=[[1,2],[2,3],[3,4],[4,5]],
                                angles=None,
                                dihedrals=[[1,2,3,4]],
                                posres=[1]
                        )
                        
C6=normal_molecule(molname='C6',atoms=[{'atomtype':'C03','resid':1,'resname':'C6','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                       {'atomtype':'C03','resid':1,'resname':'C6','q':0,'mass':calculate_mass_from_formula('C3H7')}
                                      ],
                                bonds=[[1,2]],
                                angles=None,
                                dihedrals=None,
                                posres=[1]
                      )                               

C7=normal_molecule(molname='C7',atoms=[{'atomtype':'C03','resid':1,'resname':'C7','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                       {'atomtype':'C04','resid':1,'resname':'C7','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2]],
                                angles=None,
                                dihedrals=None,
                                posres=[1]
                      ) 
                        
C8=normal_molecule(molname='C8',atoms=[{'atomtype':'C04','resid':1,'resname':'C8','q':0,'mass':calculate_mass_from_formula('C4H9')},
                                       {'atomtype':'C04','resid':1,'resname':'C8','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2]],
                                angles=None,
                                dihedrals=None,
                                posres=[1]
                      )                          

C9=normal_molecule(molname='C9',atoms=[{'atomtype':'C03','resid':1,'resname':'C9','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                       {'atomtype':'C03','resid':1,'resname':'C9','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C03','resid':1,'resname':'C9','q':0,'mass':calculate_mass_from_formula('C3H7')}
                                      ],
                                bonds=[[1,2],[2,3]],
                                angles=None,
                                dihedrals=None,
                                posres=[1]
                      )                      

C10=normal_molecule(molname='C10',atoms=[{'atomtype':'C03','resid':1,'resname':'C10','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                       {'atomtype':'C03','resid':1,'resname':'C9','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C04','resid':1,'resname':'C9','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3]],
                                angles=None,
                                dihedrals=None,
                                posres=[1]
                      )  

C11=normal_molecule(molname='C11',atoms=[{'atomtype':'C03','resid':1,'resname':'C11','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                       {'atomtype':'C04','resid':1,'resname':'C11','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C11','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3]],
                                angles=None,
                                dihedrals=None,
                                posres=[1]
                      ) 

C12=normal_molecule(molname='C12',atoms=[{'atomtype':'C04','resid':1,'resname':'C12','q':0,'mass':calculate_mass_from_formula('C4H9')},
                                       {'atomtype':'C04','resid':1,'resname':'C12','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C12','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3]],
                                angles=None,
                                dihedrals=None,
                                posres=[1]
                      )   

C13=normal_molecule(molname='C13',atoms=[{'atomtype':'C04','resid':1,'resname':'C13','q':0,'mass':calculate_mass_from_formula('C4H9')},
                                       {'atomtype':'C03','resid':1,'resname':'C13','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C03','resid':1,'resname':'C13','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C03','resid':1,'resname':'C13','q':0,'mass':calculate_mass_from_formula('C3H7')}
                                      ],
                                bonds=[[1,2],[2,3],[3,4]],
                                angles=None,
                                dihedrals=[[1,2,3,4]],
                                posres=[1]
                      )                  

C14=normal_molecule(molname='C14',atoms=[{'atomtype':'C04','resid':1,'resname':'C14','q':0,'mass':calculate_mass_from_formula('C4H9')},
                                       {'atomtype':'C04','resid':1,'resname':'C14','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C03','resid':1,'resname':'C14','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C03','resid':1,'resname':'C14','q':0,'mass':calculate_mass_from_formula('C3H7')}
                                      ],
                                bonds=[[1,2],[2,3],[3,4]],
                                angles=None,
                                dihedrals=[[1,2,3,4]],
                                posres=[1]
                      ) 

C15=normal_molecule(molname='C15',atoms=[{'atomtype':'C04','resid':1,'resname':'C15','q':0,'mass':calculate_mass_from_formula('C4H9')},
                                       {'atomtype':'C04','resid':1,'resname':'C15','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C15','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C03','resid':1,'resname':'C15','q':0,'mass':calculate_mass_from_formula('C3H7')}
                                      ],
                                bonds=[[1,2],[2,3],[3,4]],
                                angles=None,
                                dihedrals=[[1,2,3,4]],
                                posres=[1]
                      ) 

C16=normal_molecule(molname='C16',atoms=[{'atomtype':'C04','resid':1,'resname':'C16','q':0,'mass':calculate_mass_from_formula('C4H9')},
                                       {'atomtype':'C04','resid':1,'resname':'C16','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C16','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C16','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3],[3,4]],
                                angles=None,
                                dihedrals=[[1,2,3,4]],
                                posres=[1]
                      )                  
C17=normal_molecule(molname='C17',atoms=[{'atomtype':'C03','resid':1,'resname':'C17','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                       {'atomtype':'C03','resid':1,'resname':'C17','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C03','resid':1,'resname':'C17','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C04','resid':1,'resname':'C17','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C17','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3],[3,4],[4,5]],
                                angles=None,
                                dihedrals=[[1,2,3,4],[2,3,4,5]],
                                posres=[1]
                      )    

C18=normal_molecule(molname='C18',atoms=[{'atomtype':'C03','resid':1,'resname':'C18','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                       {'atomtype':'C03','resid':1,'resname':'C18','q':0,'mass':calculate_mass_from_formula('C3H6')},
                                       {'atomtype':'C04','resid':1,'resname':'C18','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C18','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C18','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3],[3,4],[4,5]],
                                angles=None,
                                dihedrals=[[1,2,3,4],[2,3,4,5]],
                                posres=[1]
                      )    

C19=normal_molecule(molname='C19',atoms=[{'atomtype':'C03','resid':1,'resname':'C19','q':0,'mass':calculate_mass_from_formula('C3H7')},
                                       {'atomtype':'C04','resid':1,'resname':'C19','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C19','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C19','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C19','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3],[3,4],[4,5]],
                                angles=None,
                                dihedrals=[[1,2,3,4],[2,3,4,5]],
                                posres=[1]
                      )    


C20=normal_molecule(molname='C20',atoms=[{'atomtype':'C04','resid':1,'resname':'C20','q':0,'mass':calculate_mass_from_formula('C4H9')},
                                       {'atomtype':'C04','resid':1,'resname':'C20','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C20','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C20','q':0,'mass':calculate_mass_from_formula('C4H8')},
                                       {'atomtype':'C04','resid':1,'resname':'C20','q':0,'mass':calculate_mass_from_formula('C4H9')}
                                      ],
                                bonds=[[1,2],[2,3],[3,4],[4,5]],
                                angles=None,
                                dihedrals=[[1,2,3,4],[2,3,4,5]],
                                posres=[1]
                      )    