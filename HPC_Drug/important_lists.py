######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contais useful lists, tuples and dictionaries
"""

#This tuple contains metal ions and some other ions that can be found in pdb structures
metals = ('AL', 'BA',
        'CA', 'CD',
        'CL', 'CO',
        'CS', 'CU',
        'CU1', 'CUA',
        'HG', 'IN',
        'IOD', 'K',
        'MG', 'MN3',
        'NA', 'PB',
        'PT', 'RB',
        'TB', 'TL',
        'WO4', 'YB',
        'ZN', 'CAL',
        'CES', 'CLA',
        'POT', 'SOD',
        'ZN2', 'FE',
        'FE2', 'MN',
        'NI', 'LI',
        'Y', 'U')

hist_resnames = ('HIS', 'HIE', 'HID', 'HSD', 'HSE')
cyst_resnames = ('CYS',' CYM', 'CYZ', 'CYSH', 'CYSM', 'CYX')

#contains the trash ligands that may be contained in the protein file
trash = ('PR', 'SO4',
        'SF4', 'TRS',
        'DMS', 'GLC',
        'PO4', 'MBO',
        'GOL', 'EDO',
        'LYS', 'GLU',
        'ACY', 'NAG'
        'ALA', 'PHE',
        'PEG', 'ACT',
        'LEU', 'FMT'
        'BU3', 'TYR',
        'CYS', 'ASN',
        'CSO', 'CSD')

# Contains the trash metal ions that need to be removed
trash_ions = ('IOD', 'K',
            'NA', 'CL',
            'LI', 'CA')


#The periodic table with atom weights dictionary
#This is part of CCP4
#http://www.ccp4.ac.uk/dist/checkout/arcimboldo/src/geometry.py
atom_weights = {
      'H'  :   1.00794,
      'He' :   4.002602,
      'Li' :   6.941,
      'Be' :   9.012182,
      'B'  :  10.811,
      'C'  :  12.0107,
      'N'  :  14.0067,
      'O'  :  15.9994,
      'F'  :  18.9984032,
      'Ne' :  20.1797,
      'Na' :  22.989770,
      'Mg' :  24.3050,
      'Al' :  26.981538,
      'Si' :  28.0855,
      'P'  :  30.973761,
      'S'  :  32.065,
      'Cl' :  35.453,
      'Ar' :  39.948,
      'K'  :  39.0983,
      'Ca' :  40.078,
      'Sc' :  44.955910,
      'Ti' :  47.867,
      'V'  :  50.9415,
      'Cr' :  51.9961,
      'Mn' :  54.938049,
      'Fe' :  55.845,
      'Co' :  58.933200,
      'Ni' :  58.6934,
      'Cu' :  63.546,
      'Zn' :  65.39,
      'Ga' :  69.723,
      'Ge' :  72.64,
      'As' :  74.92160,
      'Se' :  78.96,
      'Br' :  79.904,
      'Kr' :  83.80,
      'Rb' :  85.4678,
      'Sr' :  87.62,
      'Y'  :  88.90585,
      'Zr' :  91.224,
      'Nb' :  92.90638,
      'Mo' :  95.94,
      'Tc' :  98.0,
      'Ru' : 101.07,
      'Rh' : 102.90550,
      'Pd' : 106.42,
      'Ag' : 107.8682,
      'Cd' : 112.411,
      'In' : 114.818,
      'Sn' : 118.710,
      'Sb' : 121.760,
      'Te' : 127.60,
      'I'  : 126.90447,
      'Xe' : 131.293,
      'Cs' : 132.90545,
      'Ba' : 137.327,
      'La' : 138.9055,
      'Ce' : 140.116,
      'Pr' : 140.90765,
      'Nd' : 144.24,
      'Pm' : 145.0,
      'Sm' : 150.36,
      'Eu' : 151.964,
      'Gd' : 157.25,
      'Tb' : 158.92534,
      'Dy' : 162.50,
      'Ho' : 164.93032,
      'Er' : 167.259,
      'Tm' : 168.93421,
      'Yb' : 173.04,
      'Lu' : 174.967,
      'Hf' : 178.49,
      'Ta' : 180.9479,
      'W'  : 183.84,
      'Re' : 186.207,
      'Os' : 190.23,
      'Ir' : 192.217,
      'Pt' : 195.078,
      'Au' : 196.96655,
      'Hg' : 200.59,
      'Tl' : 204.3833,
      'Pb' : 207.2,
      'Bi' : 208.98038,
      'Po' : 208.98,
      'At' : 209.99,
      'Rn' : 222.02,
      'Fr' : 223.02,
      'Ra' : 226.03,
      'Ac' : 227.03,
      'Th' : 232.0381,
      'Pa' : 231.03588,
      'U'  : 238.02891,
      'Np' : 237.05,
      'Pu' : 244.06,
      'Am' : 243.06,
      'Cm' : 247.07,
      'Bk' : 247.07,
      'Cf' : 251.08,
      'Es' : 252.08,
      'Fm' : 257.10,
      'Md' : 258.10,
      'No' : 259.10,
      'Lr' : 262.11,
      'Rf' : 261.11,
      'Db' : 262.11,
      'Sg' : 266.12,
      'Bh' : 264.12,
      'Hs' : 269.13,
      'Mt' : 268.14,
  }


# a dictionary of how many nanoseconds per day a given processor kind
#can process on a 15000 atoms structure
#for CPU only simulations
processor_kind_ns_per_day_15000_atoms_for_cpu_only_runs = {
        'skylake' : 4.,
        'broadwell' : 5.,
        'knl' : 1.
}

#the ammount of ns per day a GPU accellerated architecture can usually generate
#on a 15000 atoms sistem
#this is a very conservative value!!!
ns_per_day_15000_on_gpu_accellerated_architectures = 15.0