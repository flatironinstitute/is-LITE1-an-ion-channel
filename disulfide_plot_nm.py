import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import seaborn as sns

traj_file = 'step7-centered.xtc'
top_file = 'step7.pdb'

traj = md.load(traj_file,top=top_file)
top = traj.topology


A_disulfide_residues = [top.select('resid %s resname CYS and name CA'%(44-35))[0],
                        top.select('resid %s resname CYS and name CA'%(149-35))[0]]
B_disulfide_residues = [top.select('resid %s resname CYS and name CA'%(44-35+388))[0],
                        top.select('resid %s resname CYS and name CA'%(149-35+388))[0]]
C_disulfide_residues = [top.select('resid %s resname CYS and name CA'%(44-35+2*388))[0],
                        top.select('resid %s resname CYS and name CA'%(149-35+2*388))[0]]
D_disulfide_residues = [top.select('resid %s resname CYS and name CA'%(44-35+3*388))[0],
                        top.select('resid %s resname CYS and name CA'%(149-35+3*388))[0]]

disulfide_A = md.compute_distances(traj,[[A_disulfide_residues[0],A_disulfide_residues[1]]])
disulfide_B = md.compute_distances(traj,[[B_disulfide_residues[0],B_disulfide_residues[1]]])
disulfide_C = md.compute_distances(traj,[[C_disulfide_residues[0],C_disulfide_residues[1]]])
disulfide_D = md.compute_distances(traj,[[D_disulfide_residues[0],D_disulfide_residues[1]]])

plt.figure(figsize=(4,2))
plt.plot(traj.time*1e-6,disulfide_A)
plt.plot(traj.time*1e-6,disulfide_B)
plt.plot(traj.time*1e-6,disulfide_C)
plt.plot(traj.time*1e-6,disulfide_D)
plt.ylabel(r'Cys44-Cys149 C$\alpha$ distance (nm)')
plt.xlabel(r'time ($\mu$s)')
plt.ylim((0,1))
sns.despine()

plt.savefig('disulfide_dist_nm.png',bbox_inches='tight',dpi=300)

combined_disulfide = np.concatenate((disulfide_A,disulfide_B,disulfide_C,disulfide_D))

rank0_disulfides = np.load('../rank0/rank0-disuflides.npy')

plt.figure(figsize=(2,2))
plt.hist(rank0_disulfides,bins=100,color='0.5');
plt.hist(combined_disulfide,bins=100,color='navy',alpha=0.8);
plt.xlabel(r'Cys44-Cys149 C$\alpha$ distance (nm)')
plt.ylabel('P')
plt.yticks([])
sns.despine()
plt.savefig('disulfide_histogram_nm.png',bbox_inches='tight',dpi=300)



