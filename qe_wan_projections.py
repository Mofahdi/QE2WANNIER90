import os
import numpy as np
import pandas as pd
from jarvis.core.kpoints import Kpoints3D
import ase
import ase.io.espresso as espres
from jarvis.core.atoms import Atoms
from ase.atoms import Atoms as AseAtoms
from jarvis.core.atoms import ase_to_atoms 


def get_struc (path, type='qe_input'):
	type=type.lower()
	if type=='poscar':
		atoms=ase.io.read(path)
	elif type=='cif':
		atoms=ase.io.read(path)
	elif type=='qe_input':
		atoms=espres.read_espresso_in(path)
	else:
		raise Exception("sorry, the code currently supports 'cif' and 'poscar' formats or 'qe_input' files. Case is ignored!")

	return atoms


def get_basis(ase_atoms): # input ase atoms
	# atomdata was obtained from the link below: (with some personal modifications)
	# https://github.com/usnistgov/jarvis/blob/master/jarvis/io/qe/outputs.py
	atomdata=dict()
	atomdata["H"] = ["s"]

	atomdata["H"] = ["s"]
	atomdata["Li"] = ["s", "p"]
	atomdata["Be"] = ["s", "p"]
	atomdata["B"] = ["s", "p"]
	atomdata["C"] = ["s", "p"]
	atomdata["N"] = ["s", "p"]
	atomdata["O"] = ["s", "p"]
	atomdata["F"] = ["s", "p"]

	atomdata["Na"] = ["s", "p"]
	atomdata["Mg"] = ["s", "p"]
	atomdata["Al"] = ["s", "p"]

	atomdata["Si"] = ["s", "p"]
	atomdata["P"] = ["s", "p"]
	atomdata["S"] = ["s", "p"]
	atomdata["Cl"] = ["s", "p"]
	atomdata["K"] = ["s", "d", "p"]
	atomdata["Ca"] = ["s", "d", "p"]
	atomdata["Sc"] = ["s", "d", "p"]
	atomdata["Ti"] = ["s", "d", "p"]
	atomdata["V"] = ["s", "d", "p"]
	atomdata["Cr"] = ["s", "d", "p"]
	atomdata["Mn"] = ["s", "d", "p"]
	atomdata["Fe"] = ["s", "d", "p"]
	atomdata["Co"] = ["s", "d", "p"]
	atomdata["Ni"] = ["s", "d", "p"]
	atomdata["Cu"] = ["s", "d", "p"]
	atomdata["Zn"] = ["s", "d", "p"]
	atomdata["Ga"] = ["s", "d", "p"]
	atomdata["Ge"] = ["s", "p"]
	atomdata["As"] = ["s", "p"]
	atomdata["Se"] = ["s", "p"]
	atomdata["Br"] = ["s", "p"]

	atomdata["Rb"] = ["s", "d", "p"]
	atomdata["Sr"] = ["s", "d", "p"]
	atomdata["Y"] = ["s", "d", "p"]
	atomdata["Zr"] = ["s", "d", "p"]
	atomdata["Nb"] = ["s", "d", "p"]
	atomdata["Mo"] = ["s", "d", "p"]
	atomdata["Tc"] = ["s", "d", "p"]
	atomdata["Ru"] = ["s", "d", "p"]
	atomdata["Rh"] = ["s", "d", "p"]
	atomdata["Pd"] = ["s", "d", "p"]
	atomdata["Ag"] = ["s", "d", "p"]
	atomdata["Cd"] = ["s", "d", "p"]
	atomdata["In"] = ["s", "d", "p"]
	atomdata["Sn"] = ["s", "p"]
	atomdata["Sb"] = ["s", "p"]
	atomdata["Te"] = ["s", "p"]
	atomdata["I"] = ["s", "p"]
	atomdata["Cs"] = ["s", "d", "p"]
	atomdata["Ba"] = ["s", "d", "p"]
	atomdata["La"] = ["s", "d"]
	atomdata["Hf"] = ["s", "d", "p"]
	atomdata["Ta"] = ["s", "d", "p"]
	atomdata["W"] = ["s", "d", "p"]
	atomdata["Re"] = ["s", "d", "p"]
	atomdata["Os"] = ["s", "d", "p"]
	atomdata["Ir"] = ["s", "d", "p"]
	atomdata["Pt"] = ["s", "d", "p"]
	atomdata["Au"] = ["s", "d", "p"]
	atomdata["Hg"] = ["s", "d", "p"]
	atomdata["Tl"] = ["s", "d", "p"]
	atomdata["Pb"] = ["s", "p"]
	atomdata["Bi"] = ["s", "p"]
	
	formula=ase_atoms.symbols.formula
	num_atoms_dict=formula.count()
	symbols=list(num_atoms_dict.keys())
	
	basis_dict=dict()
	for i, symbol in enumerate(symbols):
		basis_dict[symbol]=atomdata[symbol]
	return basis_dict, num_atoms_dict
	

def get_band_kpoints(Ase_atoms) -> list:
	atoms=ase_to_atoms(Ase_atoms, True)
	hi_kp=Kpoints3D().high_kpath(atoms=atoms)
	lines=[]; lines.append('begin kpoint_path\n')
	for i_path in range(len(hi_kp['path'])):
		for hi_symm_kpt_i, hi_symm_kpt in enumerate(hi_kp['path'][i_path]):
			i_hi_fcoord=hi_kp['kpoints'][hi_kp['path'][i_path][hi_symm_kpt_i]]
			f_hi_fcoord=hi_kp['kpoints'][hi_kp['path'][i_path][hi_symm_kpt_i+1]]
			i_hi_kpt=hi_kp['path'][i_path][hi_symm_kpt_i]
			f_hi_kpt=hi_kp['path'][i_path][hi_symm_kpt_i+1]
			#if i_hi_kpt=='\\Gamma':
			#	i_hi_kpt='G'
			i_hi_kpt= 'G' if hi_kp['path'][i_path][hi_symm_kpt_i]=='\\Gamma' else hi_kp['path'][i_path][hi_symm_kpt_i]
			f_hi_kpt= 'G' if hi_kp['path'][i_path][hi_symm_kpt_i+1]=='\\Gamma' else hi_kp['path'][i_path][hi_symm_kpt_i+1]
			hi_kpt_i=i_hi_kpt+' '+str(i_hi_fcoord[0])+' '+str(i_hi_fcoord[1])+' '+str(i_hi_fcoord[2])
			hi_kpt_f=f_hi_kpt+' '+str(f_hi_fcoord[0])+' '+str(f_hi_fcoord[1])+' '+str(f_hi_fcoord[2])
			line=hi_kpt_i+'\t'+hi_kpt_f+'\n'
			lines.append(line)
			
			if hi_symm_kpt_i==len(hi_kp['path'][i_path])-2:
				break
	lines.append('end kpoint_path\n')

	return lines
	
	
def write_wan_projections(basis_dict, num_atoms_dict, other_commands=None, ase_atoms=None,
band_points=30, projection_type='orbitals', file_name='wannier90_QE.win', SOC=True):
	f=open(file_name, mode='w')
	f.write('Begin Projections\n')
	num_wan=0
	for specie, orbitals in basis_dict.items():
		specie_proj='{specie}:'.format(specie=specie)
		specie_l_proj=specie_proj
		for orbital in orbitals:
			if orbital=='s':
				num_wan+=2*num_atoms_dict[specie]
				specie_proj=specie_proj+'s;'
				specie_l_proj=specie_l_proj+"l=0,mr=1;"
			if orbital=='p':
				num_wan+=6*num_atoms_dict[specie]
				specie_proj=specie_proj+'p;'
				specie_l_proj=specie_l_proj+"l=1,mr=1,2,3;"
			if orbital=='d':
				num_wan+=10*num_atoms_dict[specie]
				specie_proj=specie_proj+'d;'
				specie_l_proj=specie_l_proj+"l=2,mr=1,2,3,4,5;"
			if orbital=='f':
				num_wan+=14*num_atoms_dict[specie]
				specie_proj=specie_proj+'f;'
				specie_l_proj=specie_l_proj+"l=3,mr=1,2,3,4,5,6,7;"
		specie_proj=specie_proj[:-1]
		specie_l_proj=specie_l_proj[:-1]
		if projection_type=='orbitals':
			f.write('%s\n'%(specie_proj))
		elif projection_type=='qunatum_numbers':
			f.write('%s\n'%(specie_l_proj))
		else:
			raise Exception("sorry, the code currently supports 'orbitals' and 'qunatum_numbers' keywords. Case is ignored!")
	f.write('End Projections\n')

	if SOC==True:
		num_wan_used=num_wan
	elif SOC==False:
		num_wan_used=num_wan/2
	f.write('num_wann ='+str(int(num_wan_used))+'\n')

	if other_commands:
		for key, val in other_commands.items():
			f.write(key+' = '+str(val)+'\n')
			if key=='bands_plot' and val=='true':
				kpoints_path=get_band_kpoints(ase_atoms)
				f.write('bands_num_points = '+str(band_points)+'\n')
				for line in kpoints_path:
					f.write(line)
	
	f.close()


if __name__=="__main__":
	ase_atoms=get_struc('MnCo2Si_scf.in')
	basis_dict, num_atoms_dict=get_basis(ase_atoms)
	
	## you could get atoms (Jarvis atoms) then you use it as input for the "get_basis" function
	#atoms=espres.read_espresso_in('MnCo2Si_scf.in')
	#atoms=ase.io.read('POSCAR'); print(atoms)
	#atoms=ase.io.read('MnCo2Si.cif'); print(atoms.__dict__)

	wan_commands={'num_iter': 300, 'write_xyz': 'true',  'write_hr': 'true', 'guiding_centres': 'true',
	'bands_plot': 'true'}
	write_wan_projections(basis_dict, num_atoms_dict, wan_commands, ase_atoms, projection_type='orbitals', file_name='wannier90_qe.win')	