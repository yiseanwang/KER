import numpy as np
import pandas as pd
import glob
import mdtraj as md
# import matplotlib.pyplot as plt

# 先假設分裂的狀況已經決定好,所以the_frag包含fragment裡面的原子(a list) 
def kinetic_energy_final(xyz_file, the_frag, mass):
	f = open(xyz_file, 'r')  # 讀進xyz file
	lines = f.readlines()
	f.close()
	Vcx=0
	Vcy=0
	Vcz=0
	i=0
	while i < len(the_frag):
		parts = lines[int(the_frag[i])].split()  # 把第the_frag[i], 第i個原子,行分開 [x, y, z, vx, vy, vz]
		vix, viy, viz = float(parts[-3]), float(parts[-2]), float(parts[-1])
		Vcx =  Vcx + mass[i]*float(parts[-3])    # 這裡應該要有分母(m1+m2+...) 但是KER會乘回去 所以拿掉了
		Vcy =  Vcy + mass[i]*float(parts[-2])
		Vcz =  Vcz + mass[i]*float(parts[-1])
	final_K = 0.5 * 104.208*(Vcx**2+Vcy**2+Vcz**2)
    return final_K

def kinetic_energy_final(xyz_file, the_frag, masses):
	f = open(xyz_file, 'r')  # 讀進xyz file
	lines = f.readlines()
	f.close()
	M=masses.sum()
	Vcx=0
	Vcy=0
	Vcz=0
	i=0
	while i < len(the_frag):
		parts = lines[int(the_frag[i])].split()  # 把第the_frag[i], 第i個原子,行分開 [x, y, z, vx, vy, vz]
		vix, viy, viz = float(parts[-3]), float(parts[-2]), float(parts[-1])
		Vcx =  Vcx + mass[i]*(vix/M)   
		Vcy =  Vcy + mass[i]*(viy/M)
		Vcz =  Vcz + mass[i]*(viz/M)
	final_K = 0.5 * 104.208*(Vcx**2+Vcy**2+Vcz**2)
    return final_K

def velocity_of_CoM(xyz_file_path, top_file_path, masses, x, threshold=0.0):
        traj = md.load_xyz(xyz_file_path, top_file_path)#, top=None, format='xyz')
		M=masses.sum()
		Vcx =  Vcx + mass[i]*float(parts[-3])    # 這裡應該要有分母(m1+m2+...) 但是KER會乘回去 所以拿掉了
		Vcy =  Vcy + mass[i]*float(parts[-2])
		Vcz =  Vcz + mass[i]*float(parts[-1])
        distances = 10*md.compute_distances(traj, x)# [[0,1],[0,2]])
        return Vx, Vy, Vz

def potential_enegy(xyz_file_path, top_file_path, charge1, charge2):
	traj = md.load_xyz(xyz_file_path, top_file_path)
	distances = 10*md.compute_distances(traj, x)  # 要轉成bohr
	V = ((charge1 * charge2) / distances )*104.208
	return V

# Example usage:
atoms=['O1','H2','H3','O4','H5','H6']
masses = [16, 2.014, 2.014, 16, 2.014, 2.014]

files = glob.glob("*.xyz")

for ifile in files:
	xyz_file = md.load_xyz(ifile, top_file_path)
#	filename=ifile.strip(".xyz")
	datafile = open(filename+'.csv', 'w+') # func + '_' + basis + '.csv', 'w+')
	datafile.write('time')
	for i in range(1,len(masses)+1):
		datafile.write(',atom'+str(i))
	datafile.write('\n')
	kenergy_list=[]
	for i in range(0,len(masses)):
		atom_k_in_t = kinetic_energy(ifile, i)
		kenergy_list.append(atom_k_in_t)
#	print(kenergy_list)
	for i in range(0, len(kenergy_list[0])):
		datafile.write(str(i))        # time=0,1,2,3...
		for j in range(0, len(masses)):
			datafile.write(','+str(kenergy_list[j][i]))  # j for atom1, atom2,....
		datafile.write('\n')
	datafile.close()
	df=pd.read_csv(filename+'.csv')
	df['total_k']=df.iloc[:,1:7].sum(axis=1)
	df['ker']=df['total_k']-df.iloc[0,-1] # total_k 減掉一開始的k, # 應該有更好的寫法
	df.to_csv(filename+'.csv', index=False)
#print("total kinetic energy :", total_kinetic_energy)
#	print("KER :", df['ker'][-1])

mdtraj.compute_center_of_mass
