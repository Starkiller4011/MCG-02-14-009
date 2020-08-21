#!/bin/python

from numpy.linalg import svd
import numpy as n
import pylab as p
from matplotlib.mlab import PCA
import glob
import pyfits
import matplotlib
from matplotlib.ticker import *
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.stats.stats import pearsonr
import os
import subprocess
import sys


#font={'family' : 'normal','size'   : 14}
#p.rc('font', **font)

########################## Settings ####################################

#tstep=int(raw_input('Enter timestep: '))
tstep=10000
xmmbinnumber=30
nubinnumber=20
#Change between 'spectra' and 'lightcurve' for better spectral/temporal resolution
pcatype='spectra'
#Write spectra to pha files for xspec fitting
writetofile=False
#Show plot of components
plot=True
#Calculate errors (warning: increases run time by ~10x)
errorcalc=True
#Track variability in components by reconstructing spectra, produce min/max spectra for each component
reconstruct=False
minmax=True
printnorms= True
#If plotting, multiply by mean spectrum to get flux
incflux=False
ylog=False
#Apply smoothing function to spectra
smoothspec=False
#Number of spectra to plot
nspec=3
#Number of perturbed spectra to create for error calculation
nerrors=10
#If true, look for deviations between adjacent bins, rather than deviations from the mean
timescaled=False
#if true, remove files that have no counts/too short an exposure
delempty=False
#Set to True to plot variability spectrum for nspec components
varspec=False
#Plot the spectrum of the noise (only really useful for variability)
noisespec=True
bkgcorr=True
customenergies=True
energies=[0.6,1.0,1.2,1.5,1.8,2.2,2.6,3.0,3.5,4.0,4.5,5.5,6.5,8.0,10.0]
incnustar=False #### Doesn't work!
testthing=False #### No idea what this does
axisfontsize=14

emin=0.6
emax=10.0

numin=9.0
numax=40.

mult=[-1,-1,-1,1]

colours=['k','r','b','g','m','c','y']


folderlist=['30ks/']


#folderlist=['/scratch/mlparker/mrk1048/products_pca_0690870101/%ss/' % tstep,\
			#'/scratch/mlparker/mrk1048/products_pca_0690870501/%ss/' % tstep]

#folderlist=['/scratch/mlparker/mrk1048/products_pca_0150470601/%ss/' % tstep]

filetype='pha'

#filetype='fits'

if pcatype == 'lightcurve':
	filetype='lc'

if customenergies:
	nbins=len(energies)-1

if incflux and reconstruct:
	print 'Warning: reconstructed spectra are not valid if incfulx=True'
	exit()

################## Anaylsis begins here ##################################

print '\nSuzaku Principle Component Analysis, M. Parker 2014'
print 'Timestep: %ss' % tstep


rmffile = pyfits.open('spectra/fi10.rsp')
nurmf=pyfits.open('/scratch/mlparker/nustar/agn/new_mcg6/NU.rmf')
#nurmf=pyfits.open('/scratch/mlparker/nustar/agn/new_mcg6/products_3/3A_sr.rmf')

if pcatype =='spectra':
	filenumbers=[]
	specstarts=[]
	specfiles=[]
	for folder in folderlist:
		count=0
		for specfile in glob.glob('%s*fi*_g.%s' % (folder,filetype)):
			print specfile
			count+=1
			spec=pyfits.open(specfile)
			specdata=spec[1].data
			specheader=spec[1].header
			#gtidata=spec[3].data
			if 'EXPOSURE' not in specheader:
				print 'Warning! Keyword EXPOSURE not in file %s header! Deleting...' % specfile
				os.remove(specfile)
				count -=1
			else:
				#specstarts.append(gtidata[0][0])
				specfiles.append(specfile)
		filenumbers.append(count)
		if count == 0:
			print 'Warning! No files in folder %s\nIs filename correct?' % folder 
		else:
			print count,'files found'
	
	#specfiles=[s for t,s in sorted(zip(specstarts,specfiles))]
	#specstarts=sorted(specstarts)
	#for i in specfiles:
		#print i
	##print specfiles
	#exit()
	
	
	#xmm_init_time=specstarts[0]
	#specstarts=[t-xmm_init_time for t in specstarts]
	#print specstarts
	#exit()
	
	ebounds = rmffile[1].data
	
	if not customenergies:
		energies = []
		for i in n.linspace(n.log10(emin),n.log10(emax),xmmbinnumber):
			energies.append(10.**i)
	emins=[]
	emaxs=[]
	channels=[]
	for row in ebounds:
		channels.append(row[0])
		emins.append(row[1])
		emaxs.append(row[2])
	xmmbins=[]
	for e in energies:
		for i in range(1,len(channels)-1):
			if emins[i]<=e<emins[i+1]:
				xmmbins.append(channels[i])
	
	if incnustar:
		print '\nFinding simultaneous NuSTAR and XMM spectra...'
		obsids=[2,3,5]
		nuenergies=[]
		nuchans=[]
		numins=[]
		nubins=[]
		nustarts=[]
		nufiles=[]
		
		for i in n.linspace(n.log10(numin),n.log10(numax),nubinnumber):
			nuenergies.append(10.**i)
		nubounds=nurmf[1].data
		for row in nubounds:
			nuchans.append(row[0])
			numins.append(row[1])
		#print numins
		#exit()
		#print len(nuenergies)
		#print nuenergies
		for e in nuenergies:
			#print e
			for i in range(1,len(channels)-1):
				if numins[i]<=e<numins[i+1]:
					nubins.append(channels[i])
					#print e
		#print len(nubins)
		#exit()
		for specfile in glob.glob('%sobs*A_tstep%s_no*_sr.pha' % (nufolder,tstep)):
			spec=pyfits.open(specfile)
			specheader=spec[1].header
			exptime=specheader['EXPOSURE']
			if exptime >tstep/4.:
				gtidata=spec[2].data
				nustarts.append(gtidata[0][0])
				nufiles.append(specfile)
		nufiles=[s for t,s in sorted(zip(nustarts,nufiles))]
		#for nufile in nufiles:
			#print nufile
		#exit()
		
		nustarts=sorted(nustarts)
		nustar_init_time=nustarts[0]
		nustarts=[t-nustar_init_time for t in nustarts]
		
		jointstarts=[]
		for xmm_start in specstarts:
			
			for nu_start in nustarts:
				if abs(nu_start-xmm_start)<tstep/5.:
					jointstarts.append((nu_start,xmm_start))
				#else:
					#print nu_start,xmm_start,abs(nu_start-xmm_start)
		temp_nuspecs=[]
		temp_specfiles=[]
		for i in jointstarts:
			nuindex=nustarts.index(i[0])
			xmmindex=specstarts.index(i[1])
			#print nuindex,xmmindex
			temp_nuspecs.append(nufiles[nuindex])
			temp_specfiles.append(specfiles[xmmindex])
		specfiles=temp_specfiles
		nufiles=temp_nuspecs
		#print specfiles
		#print nufiles


elif pcatype=='lightcurve':
	print blah
	for folder in folderlist:
		count=0
		for lcfile in glob.glob(folder+'*%ss*lc' % tstep):
			filepars=lcfile[len(folder):].split('_')
			count+=1
			minchan=filepars[1]
			maxchan=filepars[2]
			channels.append(int(minchan))
		channels.append(int(maxchan))
		filenumbers.append(count)
	channels=sorted(channels)

else:
	print 'PCA type not recognised.'
	exit()
print 'Total intervals: ',sum(filenumbers)

############################# Functions ######################################
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'

def sumfunction(x,a,b,c):
	newvals=[]
	for xi in x:
		A=specs[0][xi]
		B=specs[1][xi]
		C=specs[2][xi]
		newvals.append(a*A+b*B+c*C)
	return newvals

def meancalc(listvals):
	newlist=[]
	for i in range(0,len(listvals)-1):
		newlist.append((listvals[i]+listvals[i+1])/2.)
	return newlist

def fitfunc(spectrum):
	x=range(0,len(spectrum))
	fit=curve_fit(sumfunction,x,spectrum)
	return fit

def meanfit():
	'''try to fit the components to the mean, to see if it can be explained purely 
	as the sum of the PCA components'''
	outspec=[]
	
	return 0.


def perturb(spectrum):
	'''Add/subtract a value between N^1/2 and 0 from each bin in a spectrum'''
	newspec=[]
	for count in spectrum:
		#randval=n.random.rand()*2.-1.
		randval=n.random.randn()
		if count>0.:
			newcount=count+randval*count**0.5
		else:
			newcount=0.
		newspec.append(max([newcount,0.]))
	return newspec

def fraction_finder(specset):
	totals=[]
	for i in range(0,len(specset[0])):
		total=0
		for j in range(0,len(specset)):
			total+=specset[j][i]
		totals.append(total)
	newspecset=[]
	for j in range(0,len(specset)):
		newspec=[]
		spec=specset[j]
		for i in range(0,len(spec)):
			fraction=spec[i]/totals[i]
			newspec.append(fraction)
		newspecset.append(newspec)
	return newspecset

def line(xs,fit):
	return [x*fit[0]+fit[1] for x in xs]

def rebin(chans,countlist,bins):
	newcounts=[]
	for i in range(0,len(bins)-1):
		newcounts.append(0)
	for i in range(0,len(chans)):
		channel=chans[i]
		count=countlist[i]
		for j in range(0,len(bins)-1):
			if bins[j]<channel<bins[j+1]:
				newcounts[j]+=count
	return newcounts

def smooth(data,order=2):
	outdata=[]
	outdata.append((data[0]+data[1])/2.)
	for i in range(1,len(data)-1):
		outdata.append((data[i-1]+data[i]+data[i+1])/3)
	outdata.append((data[-1]+data[-2])/2)
	return outdata

def percent(fractions):
	pcnts=[]
	for f in fractions:
		pcnts.append(f*100.)
	return pcnts

def writebins(binlist):
	bfname='xmm_bins_%s.txt' % len(binlist)
	print '\nWriting bins to file %s' % bfname
	if os.path.exists(bfname):
		print bfname, 'already exists. Overwiting...'
		os.remove(bfname)
	binfile=open(bfname, 'w')
	startstr= '0 '+str(min(binlist)-1)+' -1\n'
	binfile.write(startstr)
	for i in range(0,len(binlist)-1):
		line=str(binlist[i])+' '+str(binlist[i+1]-1)+' '+str(binlist[i+1]-binlist[i])+'\n'
		binfile.write(line)
	endstr=str(binlist[-1])+' '+str(channels[-1])+' '+'-1'
	binfile.write(endstr)

def toolbarupdate(fraction):
	sys.stdout.write('[')
	sys.stdout.write('%s' % ('=' * fraction))
	sys.stdout.write('%s' % (' ' * (toolbar_width-fraction)))
	sys.stdout.write(']')
	sys.stdout.write("\b" * (toolbar_width+2))
	sys.stdout.flush()
	
def readspec(specfile,first):
	
	spec=pyfits.open(specfile)
	specdata=spec[1].data
	
	if bkgcorr:
		backfile = specfile.replace('g' ,'bg')
		back=pyfits.open(backfile)
		backdata=back[1].data
	
	specheader=spec[1].header
	exptime=specheader['EXPOSURE']
	counts=[]
	fluxes=[]
	sumcounts=0.
	#if first:
	for i in range(0,len(specdata)):
		row=specdata[i]
		if first:
			xmm_channels.append(row[0])
		if bkgcorr:
			backrow=backdata[i]
			flux=row[1]-backrow[1]
		else:
			flux=row[1]
		fluxes.append(flux)
		counts.append(flux/exptime)
	if first:
		return counts,fluxes,exptime,xmm_channels
	else:
		return counts,fluxes,exptime

def steppedline(xs,ys,log=False,**kwargs):
	'''Produce a plot with a stepped line'''
	newxs = []
	newys = []
	newxs.append(xs[0])
	newys.append(ys[0])
	newxs.append(xs[1])
	newys.append(ys[0])
	for i in range(1,len(xs)-1):
		newxs.append(xs[i])
		newys.append(ys[i])
		newxs.append(xs[i+1])
		newys.append(ys[i])

	newys.append(ys[-1])
	newxs.append(xs[-1])
	if log == True:
		p.loglog(newxs,newys,**kwargs)
		#p.plot(newxs,newys,**kwargs)
	else:
		p.plot(newxs,newys,**kwargs)

def nureadspec(specfile,first):
	spec=pyfits.open(specfile)
	specdata=spec[1].data
	
	if bkgcorr:
		backfile = specfile.replace('_sr.pha','bk.pha')
		back=pyfits.open(backfile)
		backdata=back[1].data
	
	specheader=spec[1].header
	exptime=specheader['EXPOSURE']
	counts=[]
	fluxes=[]
	sumcounts=0.
	#if first:
	for i in range(0,len(specdata)):
		row=specdata[i]
		if first:
			nustar_channels.append(row[0])
		if bkgcorr:
			backrow=backdata[i]
			flux=row[1]-backrow[1]
		else:
			flux=row[1]
		fluxes.append(flux)
		counts.append(flux/exptime)
	if first:
		return counts,fluxes,exptime,nustar_channels
	else:
		return counts,fluxes,exptime
	

def readlc():
	return 0

##################################################################################
#print len(energies)-1
if writetofile:
	writebins(xmmbins)

#XMM PCA:

toolbar_width=60
rand_spectra=[]
xmm_datalist=[]
xmm_channels=[]
nustar_channels=[]
first=True
exposures=[]

if errorcalc:
	print '\nReading data and calculating perturbed spectra...\n'
else:
	print '\nReading data...\n'

foldernum=0
totalspectra=0
nfiles=len(specfiles)
#for folder in folderlist:
filecount=0
exposure=0.

totalcounts=0

for i in range(0,len(specfiles)):
	specfile=specfiles[i]
	if incnustar:
		nuspecfile=nufiles[i]
	fraction=int(toolbar_width*filecount/nfiles)
	toolbarupdate(fraction)
	
	if pcatype == 'spectra':
		if first:
			counts,fluxes,exptime,xmmchannels=readspec(specfile,first)
			if incnustar:
				nuspecfile=nufiles[i]
				nucounts,nufluxes,nuexptime,nuchannels=nureadspec(nuspecfile,first)
			
			first=False
		else:
			counts,fluxes,exptime=readspec(specfile,first)
			if incnustar:
				nucounts,nufluxes,nuexptime=nureadspec(nuspecfile,first)
		totalcounts+=sum(counts)*exptime
		counts=rebin(xmm_channels,counts,xmmbins)
		#print len(energies),len(counts),len(xmmchannels),len(xmm_channels),len(xmmbins)
		if incnustar:
			#p.plot(n.log10(nucounts))
			#p.show()
			nucounts=rebin(nustar_channels,nucounts,nubins)
			#print len(nuenergies),len(nucounts),len(nustar_channels),len(nuchannels),len(nubins)
			#p.plot(nuenergies[:-1],n.log10(nucounts))
			
			#p.show()
		#exit()
	if incnustar:
		counts=counts+nucounts
	#print len(energies),len(nuenergies)
	
	#print len(energies)
	#exit()
	#p.plot(energies[:-1]+nuenergies[:-1],n.log10(counts))
	#p.show()
	#exit()
	
	#elif pcatype == 'lightcurve':
		#print 0
	
	if errorcalc:
		rand_spec=[]
		for i in range(0,nerrors):
			spec_i=perturb(fluxes)
			spec_i=[i/exptime for i in spec_i]
			spec_i=rebin(xmm_channels,spec_i,xmmbins)
			rand_spec.append(spec_i)
	
	
	total=sum(counts)
	if total>0 and exptime>tstep/10.:
		totalspectra+=1
		xmm_datalist.append(counts)
		exposure+=(exptime)
		if errorcalc:
			rand_spectra.append(rand_spec)
	else:
		sys.stdout.write('%s' % (' '*(toolbar_width+2)))
		sys.stdout.write('%s' % ('\b'*(toolbar_width+2)))
		sys.stdout.flush()
		if total <=0:
			sys.stderr.write('Warning: No counts in time bin')
			sys.stderr.write('\n')
		elif exptime <=tstep/10.:
			sys.stderr.write('Warning: Exposure time too short')
			sys.stderr.write('\n')
		if delempty:
			sys.stderr.write('Deleting file.')
			sys.stderr.write('\n')
			os.remove(specfile)
			os.remove(backfile)
		sys.stderr.write(specfile+ '\n')
		sys.stderr.write('\n')
	filecount+=1

sys.stdout.write('[%s' % ('=' * toolbar_width))
sys.stdout.write('\n')
print 'Exposure time:', exposure
print 'Total counts:',int(totalcounts),'\n'

#total_exp_time=sum([sum(i) for i in exposures])
#print 'Total exposure time:', total_exp_time
if incnustar:
	energies=energies[:-1]+nuenergies
#exit()

#for i in xmm_datalist:
	#p.plot(n.log10(i))
#p.show()
#exit()

newlist=[]
xmm_array=n.array(xmm_datalist)
means=n.mean(xmm_array,axis=0)
meanerrors=[sigma/(2.0*(totalspectra)**0.5) for sigma in n.std(xmm_array,axis=0)]
#meanerrors=[]
#meanerrors=[m/10. for m in means]
#for m in means:
	#countstotal=m*exposure
	###print countstotal
	#meanerror=m/countstotal**0.5
	##print meanerror/m
	#meanerrors.append(meanerror)
##exit()
#ax=p.subplot(111)
#p.errorbar(energies[:-1],means,meanerrors)
#ax.set_xscale('log')
#ax.set_yscale('log')
#p.show()
#exit()

if timescaled:
	print '\nSubtracting adjacent spectra...'
	for i in range(1,len(xmm_array)):
		row1=xmm_array[i-1]
		row2=xmm_array[i]
		newrow=[2*(c2-c1)/(c1+c2) for c1,c2 in zip(row1,row2)]
		newlist.append(newrow)

else:
	print '\nSubtracting mean spectrum...'

	for row in xmm_array:
		#newrow=row-means
		newrow=[(c-m)/m for c,m in zip(row,means)]
		#newrow=[c/m for c,m in zip(row,means)]
		#p.plot(newrow)
		newlist.append(newrow)
#p.show()
#exit()
if errorcalc:
	print '\nCalculating errors...'
	pca_list=[]
	new_rand_spectra=[]
	eigenvals_ptbd=[]
	for rand_spec in rand_spectra:
		newspec=[]
		for spectrum in rand_spec:
			newspec.append([(c-m)/m for c,m in zip(spectrum,means)])
		new_rand_spectra.append(newspec)
	for i in range(0,nerrors):
		datalist=[]
		for spec_set in new_rand_spectra:
			datalist.append(spec_set[i])
		error_array=n.transpose(datalist)
		U,A,V=svd(error_array)
		U=n.transpose(U)
		if i==0:
			U0=n.copy(U)
		tempvals=[]
		for k in A:
			tempvals.append(k**2)
		eigenvals_ptbd.append(tempvals)
		for j in range(0,len(U)):
			row=U[j]
			if sum([a*b for a,b in zip(U0[j],row)])<0:
				row=[k*-1 for k in row]
				U[j]=row
		pca_list.append(U)
	errors=[]
	eigenvals_ptbd=n.array(eigenvals_ptbd)
	eigenerrs=n.std(eigenvals_ptbd,axis=0)
	#print 'Eigenvector errors:'
	#print eigenerrs
	#exit()
	#iterate over eigenvectors
	for i in range(0,len(pca_list[0])):
		errors.append([])
		mean_flux=means[i]
		
		#iterate over energies
		for e in range(0,len(pca_list[0][0])):
			vals=[]
			#iterate over spectra
			for j in range(0,nerrors):
				if incflux:
					val=pca_list[j][i][e]*mean_flux
				else:
					val=pca_list[j][i][e]
				vals.append(val)
			#print vals
			#exit()
			meanval=sum(vals)/len(vals)
			variance=0.
			for val in vals:
				variance+=(val-meanval)**2/float(nerrors-1)
			error=variance**0.5
			errors[i].append(error)

xmm_array=n.transpose(n.array(newlist))

print '\nCalculating PCA...'

U,A,V=svd(xmm_array)
U=n.transpose(U)
#NuSTAR PCA:
nu_datalist=[]
nu_channels=[]
first=True

eigenvals=[]
for i in A:
	eigenvals.append(i**2)

factor=sum(eigenvals)

if printnorms:
	print '\nEigenvalues:'
for i in range(0,len(eigenvals)):
	eigenvals[i]=eigenvals[i]/factor
	if printnorms:
		if errorcalc:
			print eigenvals[i], eigenerrs[i]/factor
		else:
			print eigenvals[i]
	
	
print '\nPercentage variability in 1st %s components:' % nspec
for i in range(0,nspec):
	print 'Eigenvector %s:' % str(i+1),str(eigenvals[i]*100)[0:6], '%'
#print 'Eigenvector 2:',str(eigenvals[1]*100)[0:6], '%'#, str(sum(eigenvals[0:2])*100)[0:6]
#print 'Eigenvector 3:',str(eigenvals[2]*100)[0:6], '%'#, str(sum(eigenvals[0:3])*100)[0:6]
print 'Remaining variability:',str(sum(eigenvals[nspec:])*100)[0:6],'%'

specs=[]
specfile=open('specfile.dat','w')
varspecs=[]
for specnum in range(0,min([50,len(U)])):#nspec):
	if specnum<len(mult):
		if smoothspec:
			spec=[(i*mult[specnum]) for i in smooth(U[specnum])]
		else:
			spec=[(i*mult[specnum]) for i in U[specnum]]
	else:
		spec=U[specnum]
	if incflux:
		spec=[F*c for F,c in zip(means,spec)]
	#print str(list(spec))[1:-1]
	specs.append(spec)
#exit()

if noisespec:
	noisespec=[]
	for i in range(0,len(U[0])):
		noisespec.append(0.)
	for specnum in range(nspec,len(eigenvals)):
		noisespec=[c0+c1 for c0,c1 in zip(noisespec,U[specnum])]
	if incflux:
		noisespec=[F*c for F,c in zip(means,noisespec)]

if varspec:
	for i in range(0,nspec):
		eigenval=eigenvals[i]
		varspecs.append([eigenval*abs(j) for j in specs[i]])
	varspecs=fraction_finder(varspecs)
	if noisespec:
		noisevarspec=[]
		for i in range(0,len(U[0])):
			noisevarspec.append(0.)
		for specnum in range(nspec,len(eigenvals)):
			eigenval=eigenvals[specnum]
			noisevarspec=[c0+abs(c1)*eigenval for c0,c1 in zip(noisevarspec,U[specnum])]

if testthing:
	normedmean=n.array([1. for i in means])
	#print normedmean
	coeffs,pcov=fitfunc(normedmean)
	#print coeffs,pcov
	#outspec=[]
	#for spectrum in specs[0:3]:
	c1=coeffs[0]
	c2=coeffs[1]
	c3=coeffs[2]
	outspec=[m-(c1*a+c2*b+c3*c) for m,a,b,c in zip(normedmean,specs[0],specs[1],specs[2])]
	p.plot(outspec)
	p.show()
	exit()

if reconstruct:
	print '\nReconstructing spectra...'
	psets=[]
	esets=[]
	
	for i in range(0,nspec):
		psets.append([])
		esets.append([])
		
	fits=[]
	pearsonrs=[]
	
	for spectrum in n.transpose(xmm_array):

		coeffs,pcov=fitfunc(spectrum)
		for i in range(0,nspec):
			psets[i].append(coeffs[i])
			esets[i].append(pcov[i][i]**0.5)
	
	print '\nLinear Fits:'
	print 'Gradient, correlation coefficient'
	
	for i in range(1,nspec):
		
		fits.append(n.polyfit(psets[0],psets[i],1))
		pearsonrs.append(pearsonr(psets[0],psets[i]))
		print fits[i-1][0],pearsonrs[i-1][0]
	
	if printnorms:
		print '\nNormalisation Values:'
		for i in range(0,len(psets[0])):
			print specstarts[i],psets[0][i],esets[0][i],psets[1][i],esets[1][i],psets[2][i],esets[2][i]
		raw_input()
	#if plot:
		#for i in range(0,nspec):
			#p.plot(specstarts,psets[i])
		#p.show()
	#exit()
	if minmax:
		minspecs=[]
		maxspecs=[]
		minerrors=[]
		maxerrors=[]
		for i in range(0,nspec):
			subplotnum=int(100*nspec+10+i+1)
			ax=p.subplot(subplotnum)
			p.errorbar(energies[:-1],means,meanerrors,color='k')
			normmin=min(psets[i])
			normmax=max(psets[i])
			specmin=[normmin*c*m+m for c,m in zip(specs[i],means)]
			minspecs.append(specmin)
			specmax=[normmax*c*m+m for c,m in zip(specs[i],means)]
			maxspecs.append(specmax)
			if errorcalc:
				comp_spec=specs[i]
				comp_error=errors[i]
				errormin=[smin*((esmean/smean)**2+((esmean/smean)**2+(ecomp/comp)**2)*(comp*smean)**2)**0.5 \
					for smin,esmean,smean,ecomp,comp \
					in zip(specmin,meanerrors,means,comp_error,comp_spec)]
				errormax=[smax*((esmean/smean)**2+((esmean/smean)**2+(ecomp/comp)**2)*(comp*smean)**2)**0.5 \
					for smax,esmean,smean,ecomp,comp \
					in zip(specmax,meanerrors,means,comp_error,comp_spec)]
				minerrors.append(errormin)
				maxerrors.append(errormax)
				if plot:
					p.errorbar(energies[:-1],specmin,errormin)
					p.errorbar(energies[:-1],specmax,errormax)
			else:
				if plot:
					p.plot(energies[:-1],specmin)
					p.plot(energies[:-1],specmax)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.xaxis.set_major_locator(FixedLocator([0.1,0.2,0.5,1,2,5,10,20,50,100]))
			ax.xaxis.set_minor_locator(FixedLocator([0.3,0.4,0.6,0.7,0.8,0.9,3,4,6,7,8,9,15,30,40,60,70,80,90]))
			ax.xaxis.set_major_formatter(ScalarFormatter())
			p.xlim(min(energies),max(energies))
			p.xlabel('Energy (KeV)')
			p.ylabel('Count Rate')
			#p.suptitle('Component %s' % str(i+1))
		if plot:
			p.show()
		


legendlist=['1st order','2nd order','3rd order','4th order','5th order','6th order']

#print len(energies)-1
if plot:
	
	print '\nPlotting eigenvectors...'
	#p.figure(0)
	#ax=p.subplot(211)
	#if varspec:
		#ax=p.subplot(221)
	#else:
		#ax=p.subplot(211)
	#print len(energies),len(specs[0])
	for specnum in range(0,nspec):
		subplotnum=100*nspec+11+specnum
		ax=p.subplot(subplotnum)
		if errorcalc:
			specfilename='comp%s_spectrum.txt' % str(specnum+1)
			specfile=open(specfilename,'w')
		if errorcalc:
			row0='E%s, Counts%s, +- \n' % (str(specnum+1),str(specnum+1))
			specfile.write(row0)
			for i in range(0,len(energies)-1):
				row=str(10.**((n.log10(energies[i])+n.log10(energies[i+1]))/2.))+', '+str(specs[specnum][i])+', '+str(errors[specnum][i])+'\n'
				specfile.write(row)
			p.errorbar(energies[:-1],specs[specnum],errors[specnum],ls='none',marker='x',color=colours[specnum])
		else:
			p.plot(energies[:-1],specs[specnum],marker='x',color=colours[specnum])
		ax.set_xscale('log')
		if ylog:
			ax.set_yscale('log')
		else:
			p.hlines(0,min(energies),max(energies),linestyles='dashed')
		ax.xaxis.set_major_locator(FixedLocator([0.1,0.2,0.5,1,2,5,10,20,50,100]))
		ax.xaxis.set_minor_locator(FixedLocator([0.3,0.4,0.6,0.7,0.8,0.9,3,4,6,7,8,9,15,30,40,60,70,80,90]))
		ax.xaxis.set_major_formatter(ScalarFormatter())
		p.xlim(min(energies),max(energies))
	p.xlabel('Energy (KeV)',fontsize=axisfontsize)
	p.ylabel('Normalised Count Rate',fontsize=axisfontsize)
	p.show()
	#exit()
	
		
	#if varspec:
		#ax2=p.subplot(222)
	#else:
		#ax2=p.subplot(212)
	ax2=p.subplot(111)
	
	p.xlabel('Eigenvector',fontsize=axisfontsize)
	p.ylabel('Fractional Variability',fontsize=axisfontsize)
	eigenvalmarker='o'
	for i in range(0,nspec):
		p.plot(i+1,[eigenvals[i]],marker=eigenvalmarker,color=colours[i],ls='-',ms=8)
	p.plot(range(nspec+1,len(eigenvals)),eigenvals[nspec:-1],marker=eigenvalmarker,ls='None',color='y',ms=5)
	#ax2.set_xscale('log')
	ax2.set_yscale('log')
	#ax2.yaxis.set_major_locator(FixedLocator([100,10,1,0.1,0.01,0.001,0.0001]))
	#ax2.yaxis.set_minor_locator(FixedLocator([50,25,5,2.5,0.5,0.25,0.05,0.025,0.005,0.0025]))
	#ax2.yaxis.set_major_formatter(ScalarFormatter())
	formatter=FuncFormatter(to_percent)
	#p.xlim(0,40)
	ax2.yaxis.set_major_formatter(formatter)
	
	p.legend(legendlist[:nspec],loc='upper right')
	p.show()
	exit()
	
	if varspec:
		ax3=p.subplot(223)
		for specnum in range(0,nspec):
			p.plot(energies[:-1],varspecs[specnum],marker='x',color=colours[specnum])
		if noisespec:
			p.plot(energies[:-1],noisevarspec,marker='x',color='y')
		ax3.set_yscale('log')
		ax3.set_xscale('log')
		ax3.xaxis.set_major_locator(FixedLocator([0.1,0.2,0.5,1,2,5,10,20,50,100]))
		ax3.xaxis.set_minor_locator(FixedLocator([0.3,0.4,0.6,0.7,0.8,0.9,3,4,6,7,8,9,15,30,40,60,70,80,90]))
		ax3.xaxis.set_major_formatter(ScalarFormatter())
		p.xlim(min(energies),max(energies))
		p.xlabel('Energy (KeV)')
		p.ylabel('Variance')
		#p.plot(range(1,10))
	
	if reconstruct:
		if varspec:
			ax4=p.subplot(224)
		else:
			ax3=p.subplot(223)
			
		xmin=min(psets[0])
		xmax=max(psets[0])
		xs=n.linspace(xmin,xmax,100)
		
		for i in range(1,nspec):
			p.errorbar(psets[0],psets[i],esets[i],esets[0],marker='x',ls='None',color=colours[i])
			if pearsonrs[i-1][0] >0.1:
				p.plot(xs,line(xs,fits[i-1]),color=colours[i])
		
		p.xlabel('Component 1 Normalisation')
		p.ylabel('Dependent Component Normalisation')
	p.show()
		



orders=range(1,nspec+1)
if writetofile:
	print '\nWriting to file...' 
	for order in orders:
		if minmax:
			minfilename='xmm_pca_comp%s_min.txt' % order
			maxfilename='xmm_pca_comp%s_max.txt' % order
			if os.path.exists(minfilename):
				print minfilename, 'already exists. Overwiting...'
				os.remove(minfilename)
			if os.path.exists(maxfilename):
				print maxfilename, 'already exists. Overwiting...'
				os.remove(maxfilename)
			minfile=open(minfilename,'w')
			maxfile=open(maxfilename,'w')
			for i in range(0,len(energies)-1):
				#print i
				elow=energies[i]
				ehigh=energies[i+1]
				emean=(ehigh+elow)/2.
				countmin=minspecs[order-1][i]
				countmax=maxspecs[order-1][i]
				minerror=minerrors[order-1][i]
				maxerror=maxerrors[order-1][i]
				minrow = str(elow)+' '+str(ehigh)+' '+str(countmin)+' '+str(minerror)+'\n'
				maxrow = str(elow)+' '+str(ehigh)+' '+str(countmax)+' '+str(maxerror)+'\n'
				minfile.write(minrow)
				maxfile.write(maxrow)
			
			subprocess.call(['./xmm_phamaker.sh',str(xmmbinnumber),'xmm_pca_comp%s_min' % order])
			subprocess.call(['./xmm_phamaker.sh',str(xmmbinnumber),'xmm_pca_comp%s_max' % order])
		else:
			outfilename='xmm_pca_comp_%s_bins_%s.txt' % (order,xmmbinnumber)
			if os.path.exists(outfilename):
				print outfilename, 'already exists. Overwiting...'
				os.remove(outfilename)
			outfile=open(outfilename,'w')
			for i in range(0,len(energies)-1):
				elow=energies[i]
				ehigh=energies[i+1]
				emean=(ehigh+elow)/2.
				count=specs[order-1][i]
				error=errors[order-1][i]
				row = str(elow)+' '+str(ehigh)+' '+str(count)+' '+str(error)+'\n'
				outfile.write(row)
			subprocess.call(['./xmm_phamaker.sh',str(xmmbinnumber),'xmm_pca_comp_%s_bins_%s' % (order,xmmbinnumber)])


print '\n'









