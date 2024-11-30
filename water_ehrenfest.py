from pyscf import gto, dft, scf
import numpy as np
import sys
import os
sys.path.append('/storage/home/hcoda1/3/ywang4107/p-jkretchmer3-0/Orlando/2024/O.P.TI.CAL./src') 
# /storage/home/hcoda1/8/vsuarez6/p-jkretchmer3-0/O.P.TI.CAL./src')
import rt_scf
import rt_utils
import basis_utils
import rt_ehrenfest
from rt_cap import MOCAP

dimer = gto.Mole()
water1 = gto.Mole()
water2 = gto.Mole()
proton = gto.Mole()

dimer.atom = '''
 O                  -1.41904607    0.10858435    0.00000006
 H                  -1.76001996   -0.36619300   -0.76029109
 H                  -1.76001996   -0.36619300    0.76029121
 O                  1.49137509   -0.00861627    0.00000100
 H                  1.86255228    0.87394752    0.00000101
 H                  0.53614403    0.12578597    0.00000098
'''
water1.atom = '''
 O                  -1.41904607    0.10858435    0.00000006
 H                  -1.76001996   -0.36619300   -0.76029109
 H                  -1.76001996   -0.36619300    0.76029121
'''
water2.atom = '''
 O                  1.49137509   -0.00861627    0.00000100
 H                  1.86255228    0.87394752    0.00000101
'''
proton.atom = '''
 H                  0.53614403    0.12578597    0.00000098
'''
water2.spin = 1
proton.spin = 1

dimer.basis = 'augccpvdz'
water1.basis = 'augccpvdz'
water2.basis = 'augccpvdz'
proton.basis = 'augccpvdz'

dimer.build()
water1.build()
water2.build()
proton.build()

# CAM-B3LYP
# dimer = dft.UKS(dimer); dimer.xc = 'CAMB3LYP' 
# water1 = dft.UKS(water1); water1.xc = 'CAMB3LYP' 
# water2 = dft.UKS(water2); water2.xc = 'CAMB3LYP' 
# proton = dft.UKS(proton); proton.xc = 'CAMB3LYP' 

# LC-PBE 
dimer  = dft.UKS(dimer);  dimer.xc  = 'HYB_GGA_XC_LC_PBEOP, -1 * GGA_C_OP_PBE + PBE';  dimer._numint.omega = 0.516;  dimer._numint.alpha = 0.0;  dimer._numint.beta = 1.0
water1 = dft.UKS(water1); water1.xc = 'HYB_GGA_XC_LC_PBEOP, -1 * GGA_C_OP_PBE + PBE'; water1._numint.omega = 0.516; water1._numint.alpha = 0.0; water1._numint.beta = 1.0
water2 = dft.UKS(water2); water2.xc = 'HYB_GGA_XC_LC_PBEOP, -1 * GGA_C_OP_PBE + PBE'; water2._numint.omega = 0.516; water2._numint.alpha = 0.0; water2._numint.beta = 1.0
proton = dft.UKS(proton); proton.xc = 'HYB_GGA_XC_LC_PBEOP, -1 * GGA_C_OP_PBE + PBE'; proton._numint.omega = 0.516; proton._numint.alpha = 0.0; proton._numint.beta = 1.0

dimer.kernel()
water1.kernel()
water2.kernel()
proton.kernel()

rt_water = rt_ehrenfest.RT_Ehrenfest(dimer, 1., 10, filename="H2O", prop="magnus_interpol", frequency=1, chkfile=None, verbose=6, Ne_step=1, N_step=1, get_mo_coeff_print = rt_utils.get_scf_orbitals)
#                                           dt,real-time, output,   , integrator,           , printing freq,           ,          , velocity rescaling,              
rt_water.nuc.mass[1]*=2
rt_water.nuc.mass[2]*=2
rt_water.nuc.mass[4]*=2
rt_water.nuc.mass[5]*=2
# Declare which observables to be calculated/printed
rt_water.observables.update(energy=True, mo_occ=True, charge=True, nuclei=True)

# Create object for complex absorbing potential and add to rt object
CAP = MOCAP(0.5, 0.0477, 1.0, 10.0)
rt_water.add_potential(CAP)

# Remove electron from 4th molecular orbital (in SCF basis)
# Input the two water fragments for their charge to be calculated
rt_utils.excite(rt_water, 4)
rt_water._scf.mo_occ = rt_water.occ
rt_utils.input_fragments(rt_water, water1, water2, proton)

rt_water.kernel()

