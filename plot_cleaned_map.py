#!/home/rmvd2/venv/bin/python3

#=============================================================================#
#                                                                             #
# NAME:     plot_cleaned_map.py.                                              #
#                                                                             #
# PURPOSE:  Plot results of Bayesian parametric foreground removal            #
#                                                                             #
# AUTHOR:   April-2021 rmvd2                                                  #
#                                                                             #
# CITE:     Please cite the relevant paper (de Belsunce et al. in prep.)      #
#=============================================================================#


import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import sys


base_dir    = '/rds/user/rmvd2/rds-stg20-cosmoshare/rmvd2/cmb/foregrounds/'
input_dir   = base_dir+'data/input/'
results_dir = base_dir+'data/results/stg_data/'
plot_dir    = base_dir+'/data/results/stg_data/plots/'

print('PLOT ALL MAPS')
print('0: I Stokes map')
print('1: Q Stokes map')
print('2: U Stokes map')

solveforamps  = False
if solveforamps:
	amps_file = '_amps_'
else:
	amps_file = ''

#num = input ('Choose which map to compute:')
for num in range(3):
	if num == 0: 
		print('I')
		stokes_map = 'I'
		min_ = -200.
		max_ = 200.
	elif num == 1:
		print('Q') 
		stokes_map = 'Q'
		min_ = -2.
		max_ =  2.
	elif num == 2: 
		print('U')
		stokes_map = 'U'
		min_ = -3.
		max_ =  3.

	amps_file = ''
	fcmbout      = results_dir+'cmbout_{}{}.dat'.format(stokes_map,amps_file)
	fgsout       = results_dir+'fgsout_{}{}.dat'.format(stokes_map,amps_file)
	fscanout     = results_dir+'amps_scan_out_{}{}.dat'.format(stokes_map,amps_file)
	fg = np.loadtxt(fgsout)
	cmb = np.loadtxt(fcmbout)
	#print(fg.shape, cmb.shape)
	amps_file = '_amps_'
	fcmbout      = results_dir+'cmbout_{}{}.dat'.format(stokes_map,amps_file)
	fgsout       = results_dir+'fgsout_{}{}.dat'.format(stokes_map,amps_file)
	fscanout     = results_dir+'amps_scan_out_{}{}.dat'.format(stokes_map,amps_file)
	fg_amps = np.loadtxt(fgsout)
	cmb_amps = np.loadtxt(fcmbout)
	#print(fg_amps.shape, cmb_amps.shape)
	min_amps_I = -200.
	max_amps_I = 200.
	min_amps_Q = -3.
	max_amps_Q = 3.
	min_amps_U = -3.
	max_amps_U = 3.

	plt.figure(1)
	hp.mollview(cmb[:,0], nest=True, min=min_, max=max_, title='cmb out {}'.format(stokes_map))
	plt.savefig(plot_dir+'cmb_out_{}.pdf'.format(stokes_map))
	plt.close(1)

	plt.figure(2)
	hp.mollview(cmb_amps[:,0], nest=True, title='cmb out {} + solve for amplitudes'.format(stokes_map))
	plt.savefig(plot_dir+'cmb_out{}{}.pdf'.format(amps_file,stokes_map))
	plt.close(2)

	plt.figure(3)
	hp.mollview(cmb[:,0]-cmb_amps[:,0], nest=True, min=min_, max=max_, title='DIFF: cmb out {} + solve for amplitudes'.format(stokes_map))
	plt.savefig(plot_dir+'cmb_out{}{}_diff.pdf'.format(amps_file,stokes_map))
	plt.close(3)

	plt.figure(4)
	hp.mollview(cmb[:,1], nest=True, title='cmb err {}'.format(stokes_map))
	plt.savefig(plot_dir+'cmb_err_{}.pdf'.format(stokes_map))
	plt.close(4)

	plt.figure(5)
	hp.mollview(cmb_amps[:,1], nest=True, title='cmb err {} + solve for amplitudes'.format(stokes_map))
	plt.savefig(plot_dir+'cmb_err{}{}.pdf'.format(amps_file,stokes_map))
	plt.close(5)

	plt.figure(6)
	hp.mollview(cmb[:,1]-cmb_amps[:,1], nest=True, title='DIFF: cmb out {} + solve for amplitudes'.format(stokes_map))
	plt.savefig(plot_dir+'cmb_err{}{}_diff.pdf'.format(amps_file,stokes_map))
	plt.close(6)

	plt.figure(7)
	plt.plot(cmb[:,2])
	plt.ylim(cmb[:10,2].min(), cmb[:10,2].max())
	plt.savefig(plot_dir+'cmb_mlnlike_{}.pdf'.format(stokes_map))
	plt.close(7)

	plt.figure(8)
	plt.plot(cmb_amps[:,2])
	plt.ylim(cmb_amps[:10,2].min(), cmb_amps[:10,2].max())
	plt.savefig(plot_dir+'cmb_mlnlike{}{}.pdf'.format(amps_file,stokes_map))
	plt.close(8)

	fg_names = np.array(['ampsync','specsync','ampdust','specdust'])
	for i in range(len(fg_names)):
		plt.figure()
		hp.mollview(fg[:,i], nest=True, title='{} - fg: {}'.format(stokes_map, fg_names[i]))
		plt.savefig(plot_dir+'fg_{}_{}.pdf'.format(fg_names[i], stokes_map))
		plt.close()
		plt.figure()
		hp.mollview(fg_amps[:,i], nest=True, title='{} - fg: {} + solve for amps'.format(stokes_map, fg_names[i]))
		plt.savefig(plot_dir+'fg_{}{}{}.pdf'.format(fg_names[i], amps_file, stokes_map))
		plt.close()

'''
print('read in scan results:')
scan_res = np.loadtxt('amps_scan_out_I.dat')
print(scan_res.shape)
bs = np.unique(scan_res[:,2])
bd = np.unique(scan_res[:,4])
action = scan_res[:,5]
X, Y = np.meshgrid(bs, bd)
Z=action.reshape(len(bd),len(bs))
plt.figure();plt.pcolormesh(X,Y,Z);plt.colorbar();plt.savefig('scan_grid.pdf');plt.show()
'''



