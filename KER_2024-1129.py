import pandas as pd
import os
import glob
import seaborn as sns
import matplotlib as plt
import numpy as np 
#------------------------------
def kinetic_energy_final(xyz_file, the_frag, masses):
	f = open(xyz_file, 'r')  # 讀進xyz file
	lines = f.readlines()
	f.close()
	M=0
	i=0
	while i < len(the_frag):
		M=M+masses[int(the_frag[i])]
		i=i+1
	Vcx=0
	Vcy=0
	Vcz=0
	i=0
	while i < len(the_frag):
		atom = lines[int(the_frag[i])].split()  # 把第the_frag[i], 第i個原子,行分開 [x, y, z, vx, vy, vz]
		vix, viy, viz = float(atom[-3]), float(atom[-2]), float(atom[-1])
		Vcx =  Vcx + mass[i]*(vix/M)   
		Vcy =  Vcy + mass[i]*(viy/M)
		Vcz =  Vcz + mass[i]*(viz/M)
	final_K = 0.5 * 104.208 * M * (Vcx**2+Vcy**2+Vcz**2)
    return final_K
#--------------------------------
def potential_enegy(xyz_file_path, top_file_path, charge1, charge2, frags):
	traj = md.load_xyz(xyz_file_path, top_file_path)
	distances = 10*md.compute_distances(traj, frags)/0.529 # 要轉成bohr
	V = ((charge1 * charge2) / distances )*104.208
	return V
#--------------------------------
def calculate_bond_distances(xyz_file_path, top_file_path, x, threshold=0.0):
	x = itertools.combinations(range(6), 2)  # all the atom pairs
	traj = md.load_xyz(xyz_file_path, top_file_path)#, top=None, format='xyz')
	distances = 10*md.compute_distances(traj, x)# [[0,1],[0,2]])
	return distances # bonds
#----------------------------------
top_file_path = "/storage/home/hcoda1/3/ywang4107/p-jkretchmer3-0/Orlando/2024/water-dimer_Joseph.mol2"
threshold = 2.5
atoms=['O1','H2','H3','O4','H5','H6']
masses = [16, 2.014, 2.014, 16, 2.014, 2.014] 
#--------------------------------
# A. 照順序吃xyz_final
# B. 吃進來的xyz_final分析鍵長
# C. 依鍵長分類fragmentation
# D. 將鍵長，分類，KER，寫進csv
# E. 拿csv來做圖
#---------------------------------
y = itertools.combinations(range(6), 2)  # all the atom pairs,這個出來是tuple??
############### 應該有更好的命名方法 ################## 
## 把tuple的逗點換成"-"
delim='-' 
x=[]
for i in y:
	print(i)
	#i="-".join(i)
	j = ''.join([str(ele) + delim for ele in i])
	j = j[ : len(j) - len(delim)]
	x.append(j)
#######################################################
csv_name = os.getcwd().split("/")[-3]
datafile = open(csv_name + ".csv", "w+")
datafile.write("filename")
############ 第一列命名 #############
for i in list(x):
	datafile.write(', '+str(i))
datafile.write(', '+ "type")
datafile.write(', '+ "KER" )
datafile.write('\n')
#----------------------------------------
# A. 照順序吃xyz_final, 生產出csv檔
files=glob.glob("*_final.xyz")
i=0
while i < len(files):
     # for ifile in files: 
	bond_len = calculate_bond_distances(files[i], top_file_path, x, threshold)
	datafile.write(files[i].strip("_final.xyz"))
	datafile.write(",")
	datafile.write(bond_len, delimiter=",",fmt='%10.3f')
    i=i+1

# B. 吃csv檔，依鍵長分類fragmentation
# C. 將鍵長，分類，KER，寫進csv
# D. 拿csv來做圖
#---------------------------------


# 1.read in the csv file
df=pd.read_csv('bond_length.csv')   #在這裡就應該加入tyep跟KER兩行

i=0
while i < len(df.index): 
    r = df.loc[i]        #r代表row的意思,這樣寫下一行比較短，或許下面改成function
## H2O+ H2O+
    if r.loc[(r[' 0-3'] > 10) & (r[' 0-1'] < 2) & (r[' 0-2'] < 2) & (r[' 3-4'] < 2)  & (r[' 3-5'] < 2)]:
        df[i,'type']= "type1"
        frag1=[0,1,2]
        frag2=[3,4,5]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+(V_1-V_0)
## H3O+ OH+
    elif r.loc[(r[' 0-3'] > 10) & (r[' 0-1'] < 2) & (r[' 0-2'] < 2) & (r[' 0-5'] < 2)  & (r[' 3-4'] < 2)]:
        df[i,'type']= "type2-a"
        frag1=[0,1,2,5]
        frag2=[3,4]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+(V_1-V_0)
    elif r.loc[(r[' 0-3'] > 10) & (r[' 0-1'] < 2) & (r[' 0-2'] < 2) & (r[' 0-4'] < 2)  & (r[' 3-5'] < 2)]:
        df[i,'type']= "type2-b"
        frag1=[0,1,2,4]
        frag2=[3,5]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+(V_1-V_0)
    elif r.loc[(r[' 0-3'] > 10) & (r[' 0-1'] < 2) & (r[' 2-3'] < 2) & (r[' 3-4'] < 2)  & (r[' 3-5'] < 2)]:
        df[i,'type']= "type2-c"
        frag1=[0,1]
        frag2=[3,2,4,5]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+(V_1-V_0)
    elif r.loc[(r[' 0-3'] > 10) & (r[' 0-2'] < 2) & (r[' 1-3'] < 2) & (r[' 3-4'] < 2)  & (r[' 3-5'] < 2)]:
        df[i,'type']= "type2-d"
        frag1=[0,2]
        frag2=[3,1,4,5]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+(V_1-V_0)
## H2O+ OH+ H | H2O+ OH H+
#                 #O1-O4              #O1-H2           #O1-H3             #O4-H5             #O4-H6          # O1-H6
    elif r.loc[(r[' 0-3'] > 10) & (r[' 0-1'] < 2) & (r[' 0-2'] < 2) & (r[' 3-4'] < 2)  & (r[' 3-5'] >3) & (r[' 0-5'] >3)]:
        df[i,'type']= "type3-a"
        frag1=[0,1,2]
        frag2=[3,4]
        frag3=[5]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        KER_3=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag3, masses)
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+KER_3+(V_1-V_0)
    elif r.loc[(r[' 0-3'] > 10) & (r[' 0-1'] < 2) & (r[' 0-2'] < 2) & (r[' 3-5'] < 2)  & (r[' 3-4'] >3) & (r[' 0-4'] >3)]:
        df[i,'type']= "type3-b"
        frag1=[0,1,2]
        frag2=[3,5]
        frag3=[4]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        KER_3=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag3, masses)
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+KER_3+(V_1-V_0)
#                 #O1-O4              #O4-H5            #O4-H6             #O1-H2          #O1-H3          # O1-H6
    elif r.loc[(r[' 0-3'] > 10) & (r[' 3-4'] < 2) & (r[' 3-5'] < 2) & (r[' 0-1'] < 2)  & (r[' 0-2'] >3) & (r[' 2-3'] >3)]:
        df[i,'type']= "type3-c"
        frag1=[0,1]
        frag2=[3,4,5]
        frag3=[2]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        KER_3=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag3, masses)
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+KER_3+(V_1-V_0)
    elif r.loc[(r[' 0-3'] > 10) & (r[' 3-4'] < 2) & (r[' 3-5'] < 2) & (r[' 0-2'] < 2)  & (r[' 0-1'] >3) & (r[' 1-3'] >3)]:
        df[i,'type']= "type3-d"
        frag1=[0,2]
        frag2=[3,4,5]
        frag3=[1]
        KER_1=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag1, masses)       # BOMD最後 frag1 
        KER_2=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag2, masses)       # BOMD最後 frag2
        KER_3=kinetic_energy_final(str(df.loc[i,0]+"_final.xyz"), frag3, masses)
        OandO=[0,3]
        V_1=potential_enegy(str(df.loc[i,0]+"_init.xyz"), top_file_path, charge1, charge2, OandO)  # Ehrenfest最後 
        V_0=potential_enegy(str(df.loc[i,0]+"_Eh0fs.xyz"), top_file_path, charge1, charge2, OandO) # Ehrenfest開始
        df[i,'KER']=KER_1+KER_2+KER_3+(V_1-V_0)
    else:
        df[i,'type']= "NAN"
# version 2
#    df['type'].apply(lambda x:'type1' if df[' 0-1'] < 2) & (df[' 0-2'] < 2) & (df[' 3-4'] < 2)  & (df[' 3-5'] < 2)]
#    else ('type2' if df[' 0-1'] < 2) & (df[' 0-2'] < 2) & (df[' 3-4'] < 2)  & (df[' 3-5'] < 2))
#    else ('type3' if df[' 0-1'] < 2) & (df[' 0-2'] < 2) & (df[' 3-4'] < 2)  & (df[' 3-5'] < 2))



print(os.getcwd().split("/")[-3])
print("H2O+ H2O +: ",count1)
print("H3O+ OH+ :", count2)
print("H2O+ OH+ H : ", count3)
#cond3=cond2.describe() #3. count it 計算
#cond3.to_csv('count.csv') # index=False可以移除第一行

# pandas自帶plot功能 (line, bar scatter, histogram, boxplot, area, pie chart)
df.plot.hist(bins=20)

