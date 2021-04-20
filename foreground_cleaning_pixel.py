#!/home/rmvd2/venv/bin/python3

#=============================================================================#
#                                                                             #
# NAME:     foreground_cleaning_pixel.py                                      #
#                                                                             #
# PURPOSE:  Perform foreground cleaning using Bayesian parametric fitting     #
#                                                                             #
# AUTHOR:   April-2021 rmvd2                                                  #
#                                                                             #
# CITE:     Please cite the relevant paper (de Belsunce et al. in prep.)      #
#=============================================================================#

import numpy as np
import time, sys

############################################################
# PRINT SETTING OF CODE
############################################################
print('==================================', flush=True)
print('   FOREGROUND CLEANING SETTINGS   ', flush=True)
print('==================================', flush=True)
print('pixel per pixel cleaning')
print('use STG data')
print('0: I Stokes map')
print('1: Q Stokes map')
print('2: U Stokes map')
num = int(sys.argv[1])
print('Map chosen: ', num)
print('==================================', flush=True)

############################################################
# PARAMS
############################################################
NSIDE = 16
NPIX  = 12*NSIDE**2
NFREQ = 7
NVAR  = 5
NAMPS = 3

# choose mode to run
solveforampsfirst  = False
solveforampstoinit = False
verbose = False
do_scan = False

#spectral index synchrotron
b_s_min  = -3.5
b_s_max  = -2.5
b_s_step = 0.02
syncspecmean = -3.0
b_s_std  =  0.3
invsyncspecvar  =  1./(b_s_std*b_s_std)

#spectral index dust
b_d_min  = 1.0
b_d_max  = 2.0
b_d_step = 0.001
dustspecmean = 1.5
b_s_std  =  0.5
invdustspecvar  =  1./(b_s_std*b_s_std)

freqs  = np.array([30.,44.,70.,100.,143.,217.,353.])
conv   = np.array([1.02348093053,1.05105460096,1.13324998912, 1.28673961626,1.65391841275,2.99242038479,12.9154764837])
e      = np.ones(NFREQ)
freq_sync = 30.
freq_dust = 353.

#constants
G_fac  = 10.**9
T_cmb  = 2.75 # Kelvin
h      = 4.135667696E-15
k      = 8.617333262E-5
z_deg_rad = 2*np.pi/360. 
thsmooth  = 5.*z_deg_rad*0.425
'''
# to compute conversion factors
A    = 100.        # muK (amplitude)
beta = 1.5         # er (spectral index)
T    = 18.         # K (T_d)
nu_0 = 857.        # GHz (reference frequency)
h    = 6.62607e-34 # Planck's constant
k_b  = 1.38065e-23 # Boltzmanns constant
Tcmb = 2.7255      # K CMB Temperature
'''

# compare cmbold and new computation
tol = 1.e-3
initcmb = 3 # 70GHz as initialisation CMB value
nsteps  = 200 # iterations per pixel
errcount= 0
converged_pixels= 0

############################################################
# SETTINGS
############################################################
if do_scan:
	scan_file = '_scan_'
else:
	scan_file=''
if solveforampsfirst and solveforampstoinit:
	amps_file = '_amps_'
else:
	amps_file = ''
if num == 0: 
	print('I')
	stokes_col_noise = 0
	stokes_col_map   = 0
	stokes_map = 'I'
elif num == 1:
	print('Q') 
	stokes_col_noise = 1
	stokes_col_map   = 1
	stokes_map = 'Q'
elif num == 2: 
	print('U')
	stokes_col_noise = 1
	stokes_col_map   = 2
	stokes_map = 'U'


############################################################
# INPUT AND OUTPUT DATA
############################################################
# maps
code_dir ='/home/rmvd2/cmb/foreground/codes/'
print('code based in: ', code_dir)
base_dir = '/rds/user/rmvd2/rds-stg20-cosmoshare/rmvd2/cmb/foregrounds/'
input_dir = base_dir+'data/input/'
results_dir = base_dir+'data/results/stg_data/'

f030=input_dir+'coaddedmap_n30Ghz.dat'
f044=input_dir+'coaddedmap_n44Ghz.dat'
f070=input_dir+'coaddedmap_n70Ghz.dat'
f100=input_dir+'coaddedmap_n100Ghz.dat'
f143=input_dir+'coaddedmap_n143Ghz.dat'
f217=input_dir+'coaddedmap_n217Ghz.dat'
f353=input_dir+'coaddedmap_n353Ghz.dat'
fdata=np.array([f030,f044,f070,f100,f143,f217,f353])

fnoise=input_dir+'noiselevels.dat'

fcmbout      = results_dir+'cmbout_{}{}{}.dat'.format(stokes_map,scan_file,amps_file)
fgsout       = results_dir+'fgsout_{}{}{}.dat'.format(stokes_map,scan_file,amps_file)
fsyncampout  = results_dir+'syncampout_{}{}{}.dat'.format(stokes_map,scan_file,amps_file)
fsyncpsout   = results_dir+'syncspout_{}{}{}.dat'.format(stokes_map,scan_file,amps_file)
fdustampout  = results_dir+'dustampout_{}{}{}.dat'.format(stokes_map,scan_file,amps_file)
fdustspout   = results_dir+'dustspout_{}{}{}.dat'.format(stokes_map,scan_file,amps_file)
fnoisematout = results_dir+'noisemat_{}{}{}.bin'.format(stokes_map,scan_file,amps_file)
fscanout     = results_dir+'amps_scan_out_{}{}{}.dat'.format(stokes_map,scan_file,amps_file)

############################################################
# initialize arrays
############################################################
amps_idx = np.array([0,1,3]) # indices for amplitudes
spec_idx = np.array([2,4]) # indices for spectral index

x        = np.zeros(NVAR)
dx       = np.zeros(NVAR)
model    = np.zeros(NVAR)
ninvdiff = np.zeros(NFREQ)
ninve    = np.zeros(NFREQ)

cmb      = np.zeros(NPIX)
cmbnoise = np.zeros(NPIX)
mlnlike  = np.zeros(NPIX)
ampsync  = np.zeros(NPIX)
specsync = np.zeros(NPIX)
ampdust  = np.zeros(NPIX)
specdust = np.zeros(NPIX)
fullnoise= np.zeros((NPIX,NVAR,NVAR))
noisemat = np.zeros((NFREQ,NFREQ))
dat      = np.zeros((NFREQ,NPIX))
# dsyncda:  snychrotron amplitude derivative
# dsyncdsp: snychrotron spectral index derivative
# ddustda:  dust amplitude derivative
# ddustdsp: dust spectral index derivative

############################################################
# DEFINE FUNCTIONS FOR COMPUTATIONS
############################################################
def s(nu):
	s = A*(nu/nu_0)**(beta+1)*((np.exp(h*nu_0/k_b/T)-1.)/(np.exp(h*nu/k_b/T)-1.))
	return s

def makefirstderiv(x, ninvdiff, dsyncda, dsyncdsp, ddustda, ddustdsp):
	fd = np.zeros(NVAR)
	fd[0] = -np.dot(e,        ninvdiff)
	fd[1] = -np.dot(dsyncda,  ninvdiff)
	fd[2] = -np.dot(dsyncdsp, ninvdiff)+invsyncspecvar*(x[2]-syncspecmean)
	fd[3] = -np.dot(ddustda,  ninvdiff)
	fd[4] = -np.dot(ddustdsp, ninvdiff)+invdustspecvar*(x[4]-dustspecmean)
	return fd

def makesecondderiv(x, etninve, ninvdiff, ninve, ninv, dsyncda, dsyncdsp, ddustda, ddustdsp, d2syncdadsp, d2syncdsp2, d2dustdadsp, d2dustdsp2):
	sd = np.zeros((NVAR, NVAR))

	sd[0,0]= etninve
	sd[0,1]= np.dot(dsyncda,  ninve)
	sd[0,2]= np.dot(dsyncdsp, ninve)
	sd[0,3]= np.dot(ddustda,  ninve)
	sd[0,4]= np.dot(ddustdsp, ninve)

	sd[1,1]= np.dot(dsyncda, np.dot(ninv,dsyncda))
	sd[1,2]=-np.dot(d2syncdadsp,ninvdiff)+ np.dot(dsyncda, np.dot(ninv,dsyncdsp))
	sd[1,3]= np.dot(dsyncda, np.dot(ninv,ddustda))
	sd[1,4]= np.dot(dsyncda, np.dot(ninv,ddustdsp))

	sd[2,2]=-np.dot(d2syncdsp2,ninvdiff)+invsyncspecvar + np.dot(dsyncdsp, np.dot(ninv,dsyncdsp))
	sd[2,3]= np.dot(dsyncdsp, np.dot(ninv,ddustda))
	sd[2,4]= np.dot(dsyncdsp, np.dot(ninv,ddustdsp))
	
	sd[3,3]= np.dot(ddustda, np.dot(ninv,ddustda))
	sd[3,4]=-np.dot(d2dustdadsp,ninvdiff) + np.dot(ddustda, np.dot(ninv,ddustdsp))

	sd[4,4]=-np.dot(d2dustdsp2,ninvdiff)+invdustspecvar+np.dot(ddustdsp, np.dot(ninv,ddustdsp))

	# mirror matrix --> MUST BE LOWER TRIANGULAR MATRIX == 0. 
	sd = sd + sd.T - np.diag(np.diag(sd))
	return sd

def makeaction(x,pixdat, model, ninvdiff):
	return 0.5*(np.dot((pixdat-model), ninvdiff)+(x[2]-syncspecmean)*(x[2]-syncspecmean)*invsyncspecvar+(x[4]-dustspecmean)*(x[4]-dustspecmean)*invdustspecvar)

def makemodel(x):
	return x[0]+x[1]*(freqs/freq_sync)**x[2]*(conv/conv[0])+x[3]*(freqs/freq_dust)**x[4]*(conv/conv[NFREQ-1])

def makeninvdiff(ninv, pixdat, model):
	return np.dot(ninv, (pixdat-model))

def makefgdiffs(x):
	dsyncda    =(freqs/freq_sync)**x[2]*conv/conv[0]
	dsyncdsp   =np.log(freqs/freq_sync)*x[1]*dsyncda
	ddustda    =(freqs/freq_dust)**x[4]*conv/conv[NFREQ-1]
	ddustdsp   =np.log(freqs/freq_dust)*x[3]*ddustda
	d2syncdadsp=np.log(freqs/freq_sync)*dsyncda
	d2syncdsp2 =np.log(freqs/freq_sync)*dsyncdsp
	d2dustdadsp=np.log(freqs/freq_dust)*ddustda
	d2dustdsp2 =np.log(freqs/freq_dust)*ddustdsp
	return dsyncda, dsyncdsp, ddustda, ddustdsp, d2syncdadsp, d2syncdsp2, d2dustdadsp, d2dustdsp2

def initx(x, pixdat, initcmb, ninv):
	#print('initialize with freq = {}'.format(freqs[initcmb]))
	x[0]=(pixdat[initcmb]) #CMB = 70 GHz
	x[1]=(pixdat[0]-pixdat[initcmb]) #sync = 30 GHz - CMB
	x[2]=syncspecmean
	x[3]=(pixdat[NFREQ-1]-pixdat[initcmb]) #dust = 353 GHz - CMB
	x[4]=dustspecmean
	if solveforampstoinit:
		model       = makemodel(x)
		ninvdiff    = makeninvdiff(ninv, pixdat, model)
		dsyncda, dsyncdsp, ddustda, ddustdsp, d2syncdadsp, d2syncdsp2, d2dustdadsp, d2dustdsp2 = makefgdiffs(x)
		x[amps_idx] = updateampsofx(dsyncda, ddustda, ninv, pixdat)
	return x

def makea(dsyncda, ddustda):
	a = np.zeros((NFREQ, NAMPS))
	a[:,0] = 1.
	a[:,1] = dsyncda
	a[:,2] = ddustda
	return a

def updateampsofx(dsyncda, ddustda, ninv, pixdat):
	# computes signal using analytic solution for singal 
	# s = (A.T N**-1 A)**-1 A.T N**-1 d
	a            = makea(dsyncda, ddustda) # NFREQxNAMPS
	atninv       = np.dot(a.T, ninv) # NAMPSxNFREQ , NFREQxNFREQ
	atninva      = np.dot(atninv,a)  # NAMPSxNFREQ , NFREQxNAMPS
	atninvainv   = np.linalg.inv(atninva) # NAMPSxNAMPS
	atninvpixdat = np.dot(atninv, pixdat) # NAMPSxNFREQ , NFREQx1
	amps         = np.dot(atninvainv, atninvpixdat) # NAMPSxNAMPS , NAMPSx1
	return amps # only first three components (amplitudes) of x - NOT spec indices

def makedx(x):
	# compute displacement for Newton-Raphson with: dx = -Hessian**-1 * gradient
	model    = makemodel(x)
	ninvdiff = makeninvdiff(ninv, pixdat, model)
	action   = makeaction(x,pixdat, model, ninvdiff)
	dsyncda, dsyncdsp, ddustda, ddustdsp, d2syncdadsp, d2syncdsp2, d2dustdadsp, d2dustdsp2 = makefgdiffs(x)
	if verbose: print('Before solving for amplitudes, x: ', x)
	if solveforampsfirst:
		x[amps_idx] = updateampsofx(dsyncda, ddustda, ninv, pixdat)
		model       = makemodel(x)
		ninvdiff    = makeninvdiff(ninv, pixdat, model)
		action      = makeaction(x,pixdat, model, ninvdiff)
		dsyncda, dsyncdsp, ddustda, ddustdsp, d2syncdadsp, d2syncdsp2, d2dustdadsp, d2dustdsp2 = makefgdiffs(x)
		if verbose: print('Before solving for amplitudes, x: ', x)
		if verbose: print('Before solving for amplitudes, model: ', model)
	fd = makefirstderiv(x, ninvdiff, dsyncda, dsyncdsp, ddustda, ddustdsp)
	sd = makesecondderiv(x, etninve, ninvdiff, ninve, ninv, dsyncda, dsyncdsp, ddustda, ddustdsp, d2syncdadsp, d2syncdsp2, d2dustdadsp, d2dustdsp2)
	sdinv = np.linalg.inv(sd)
	dx    = np.dot(-sdinv, fd)
	#print('model:', model)
	#print('ninvdiff:', ninvdiff)
	#print('action:', action)
	#print('dsyncda:', dsyncda)
	#print('dsyncdsp:', dsyncdsp)
	#print('ddustda:', ddustda)
	#print('ddustdsp:', ddustdsp)
	#print('d2syncdadsp:', d2syncdadsp)
	#print('d2syncdsp2:', d2syncdsp2)
	#print('d2dustdadsp:', d2dustdadsp)
	#print('d2dustdsp2:', d2dustdsp2)
	#print('fd:', fd)
	#print('sd:', sd)
	#print('sdinv:', sdinv)
	#print('dx:', dx)
	return dx, sdinv, action


############################################################
# START COMPUTATION
############################################################
print('read noiselevels')
noise = np.loadtxt(fnoise, usecols=stokes_col_noise)
print(noise)

noisemat += np.diag(noise*noise)
ninv = np.linalg.inv(noisemat)
ninve = np.dot(ninv, e)
etninve = np.dot(e, ninve)
print('etninve={}, corresponding rms noise={}.'.format(etninve,1./np.sqrt(etninve)))
print('noise: {}'.format(noise))

print('read maps:')
print(fdata)
for map_idx, map_names in enumerate(fdata):
	dat[map_idx, :] = np.loadtxt(map_names, usecols=stokes_col_map)


# perform calculation
print('Start calculation:')
for pixnum in range(NPIX):
	if pixnum % 500 == 0: print('pixel done: {}'.format(pixnum))
	pixdat = dat[:,pixnum]
	#print('data:', pixdat)
	x = initx(x, pixdat, initcmb, ninv)
	dx = np.zeros(NVAR)
	cmbold = x[0]
	for i in range(nsteps):
		x += dx
		#print('x:', x)
		dx, sdinv, action = makedx(x)
		if (i>5) and (x[0]>500. or x[0]<-500.):
			print('ERROR: problem with pixel {} in step {}'.format(pixnum, i))
			cmb[pixnum]=0.
			cmbnoise[pixnum]=np.sqrt(1./etninve)
			mlnlike[pixnum]=-1.
			errcount+=1
			break
		if (i>5) and ((x[0]-cmbold)<tol):
			# min five iterations per pixel
			converged_pixels+=1
			break
		cmbold=x[0]
	if (sdinv[0,0]<0.):
		print('ERROR: problem with sdinv in pixel {}'.format(pixnum))
		cmb[pixnum]=0.
		cmbnoise[pixnum]=np.sqrt(1./etninve)
		mlnlike[pixnum]=-1.
		errcount+=1
		continue
	cmb[pixnum]     = x[0]
	cmbnoise[pixnum]= np.sqrt(sdinv[0][0])
	mlnlike[pixnum] = action
	ampsync[pixnum] = x[1]
	specsync[pixnum]= x[2]
	ampdust[pixnum] = x[3]
	specdust[pixnum]= x[4]
	fullnoise[pixnum,:,:]=sdinv[:,:]

print('number of faulty pixels:    ', errcount)
print('number of converged pixels: ', converged_pixels)

print('Save outputs to:')
print(fcmbout)
print(fgsout)
print(fnoisematout)

np.savetxt(fcmbout, np.c_[cmb, cmbnoise, mlnlike])
np.savetxt(fgsout, np.c_[ampsync,specsync,ampdust,specdust])
fullnoise.tofile(fnoisematout)

if do_scan:
	print('now run on grid:')
	print('save scan results to: ', fscanout)

	scan_results = []
	for spec_s in np.arange(b_s_min, b_s_max, b_s_step):
		for spec_d in np.arange(b_d_min, b_d_max, b_d_step):
			x[spec_idx] = spec_s, spec_d
			dsyncda, dsyncdsp, ddustda, ddustdsp, d2syncdadsp, d2syncdsp2, d2dustdadsp, d2dustdsp2 = makefgdiffs(x)
			if verbose: print('Before solving for amplitudes, x: ', x)
			x[amps_idx] = updateampsofx(dsyncda, ddustda, ninv, pixdat)
			model       = makemodel(x)
			ninvdiff    = makeninvdiff(ninv, pixdat, model)
			action      = makeaction(x,pixdat, model, ninvdiff)
			if verbose: print('After solving for amplitudes, x: ', x, action)
			scan_results.append(np.concatenate((x,np.array([action]))))

	scan_results = np.stack(scan_results, axis=0)
	np.savetxt(fscanout, scan_results)

print('===============DONE===============', flush=True)


