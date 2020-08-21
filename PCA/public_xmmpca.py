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

currentdir=str(os.getcwd())+'/'

########################## Settings ####################################

tstep=int(raw_input('Enter timestep: '))
xmmbinnumber=int(raw_input('Enter desired number of bins: '))
#Write spectra to pha files for xspec fitting (Not sure if this works)
writetofile=False
#Show plot of components
plot=True
#Calculate errors (warning: increases run time by ~10x)
#This was more of an issue before I optimised it. Go nuts!
errorcalc=True
#Track variability in components by reconstructing spectra, produce min/max spectra for each component (seems to be working)
reconstruct=False
minmax=False
printnorms= False
#If plotting, multiply by mean spectrum to get flux - Probably don't do this
incflux=False
ylog=False #Why would you?
#Apply smoothing function to spectra. Again, don't bother.
smoothspec=False
#Number of spectra to plot
nspec=int(raw_input('How many spectra should be plotted? '))
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
# subtract the background spectrum. IMPORTANT: This currently assumes that the background and source regions are the same size. Feel free to fix it yourself. Also background subtraction may sometimes generate other components for some reason.
bkgcorr=False
#Change the values in the energies list and set customenergies=True to use your own binning.
customenergies=False
energies=[0.4,1.0,3.0,5.0,7.0,10.0]
axisfontsize=14
folder_stem='products_pca' #Tells the script which folders (within working directory) to look for files in
# This version assumes that the folders contain subfolders for each timestep of interest
file_stem='src'
# The RMF file is needed to convert between energy bins/channels. Any rmf for the same instrument will do, the response itself isn't used. You may have to change some other stuff to make it work with other instruments, not all fits files are created equal. I forget what to change though...
rmffile = pyfits.open('/scratch/mlparker/nustar/agn/new_mcg6/xmm/PN.rmf')

filetype='fits'

#Energy limits
emin=0.4
emax=10.0

#This list multiplies the components by +/-1. The y-axis is arbitrary, so it may be desirable to invert the components.
mult=[-1,-1,-1,1]

colours=['k','r','b','g','m','c','y']

# Find relevant folders
folderlist=[currentdir+f+'/%ss/' % str(tstep) for f in glob.glob(folder_stem+'*')]


if customenergies:
	nbins=len(energies)-1

if incflux and reconstruct:
	print 'Warning: reconstructed spectra are not valid if incfulx=True'
	exit()

################## Anaylsis begins here ##################################

print '\nXMM Principle Component Analysis, M. Parker 23/06/14'
print 'Timestep: %ss' % tstep

# RMF file is used to find relation between channel/energy

filenumbers=[]
specstarts=[]
specfiles=[]
for folder in folderlist:
	count=0
	#Locate and sort all files in given folder
	for specfile in glob.glob(folder+'*'+file_stem+'*'+filetype):
		#print specfile
		count+=1
		spec=pyfits.open(specfile)
		specdata=spec[1].data
		specheader=spec[1].header
		gtidata=spec[3].data
		if 'EXPOSURE' not in specheader:
			print 'Warning! Keyword EXPOSURE not in file %s header!' 
			if delempty:
				print 'Deleting...' % specfile
				os.remove(specfile)
			count -=1
		else:
			specstarts.append(gtidata[0][0])
			specfiles.append(specfile)
	filenumbers.append(count)
	if count == 0:
		print 'Warning! No files in folder %s\nIs filename correct?' % folder 
	#else:
		#print count,'files found'

specfiles=[s for t,s in sorted(zip(specstarts,specfiles))]
specstarts=sorted(specstarts)

xmm_init_time=specstarts[0]
specstarts=[t-xmm_init_time for t in specstarts]

ebounds = rmffile[2].data

# Find logarithmic energy bins
if not customenergies:
	energies = []
	for i in n.linspace(n.log10(emin),n.log10(emax),xmmbinnumber):
		energies.append(10.**i)
emins=[]
emaxs=[]
channels=[]
# Find channels corresponding to energy bins
for row in ebounds:
	channels.append(row[0])
	emins.append(row[1])
	emaxs.append(row[2])
xmmbins=[]
for e in energies:
	for i in range(1,len(channels)-1):
		if emins[i]<=e<emins[i+1]:
			xmmbins.append(channels[i])
	

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
	# Simple fuction for fitting a sum of PCs
	newvals=[]
	for xi in x:
		A=specs[0][xi]
		B=specs[1][xi]
		C=specs[2][xi]
		newvals.append(a*A+b*B+c*C)
	return newvals

def meancalc(listvals):
	'''Return list of mean values - I think this is used for plotting'''
	newlist=[]
	for i in range(0,len(listvals)-1):
		newlist.append((listvals[i]+listvals[i+1])/2.)
	return newlist

def fitfunc(spectrum):
	x=range(0,len(spectrum))
	fit=curve_fit(sumfunction,x,spectrum)
	return fit

def perturb(spectrum):
	'''Add/subtract a value between N^1/2 and 0 from each bin in a spectrum'''
	newspec=[c+n.random.randn()*c**0.5 for c in spectrum]
	return newspec

def fraction_finder(specset):
	'''I'm not sure what this does'''
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
	'''This is a straight line. Not sure why.'''
	return [x*fit[0]+fit[1] for x in xs]

def rebin(chans,countlist,bins):
	'''Rebinning function for spectra'''
	newcounts=[]
	newcounts=[0]*(len(bins)-1)
	binnum=0
	for i in range(min(bins),max(bins)):
		channel=chans[i]
		count=countlist[i]
		while bins[binnum+1]<=channel:
			binnum+=1
		newcounts[binnum]+=count
	return newcounts

def smooth(data,order=2):
	'''Smoothing function. I wouldn't bother, just bin more.'''
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
	'''Looks like this writes bins. Because reasons. I can figure this out, email me if you need it. Don't know why you would, when nobody knows what it does, but still...'''
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
	'''Self explanatory. Why did I bother typing that?'''
	sys.stdout.write('[')
	sys.stdout.write('%s' % ('=' * fraction))
	sys.stdout.write('%s' % (' ' * (toolbar_width-fraction)))
	sys.stdout.write(']')
	sys.stdout.write("\b" * (toolbar_width+2))
	sys.stdout.flush()
	
def readspec(specfile,first):
	'''Read spectrum from file'''
	spec=pyfits.open(specfile)
	specdata=spec[1].data
	
	if bkgcorr:
		backfile = specfile.replace('src' ,'bkg')
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

##################################################################################
if writetofile:
	writebins(xmmbins) # I think this is for making the PCs xspec friendly

toolbar_width=60
rand_spectra=[]
xmm_datalist=[]
xmm_channels=[]
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
	# Read data from the spectra files
	specfile=specfiles[i]
	
	fraction=int(toolbar_width*filecount/nfiles)
	toolbarupdate(fraction)
	
	if first:
		#If first, get channels as well as spectrum
		counts,fluxes,exptime,xmmchannels=readspec(specfile,first)
		first=False
	else:
		counts,fluxes,exptime=readspec(specfile,first)
	
	totalcounts+=sum(counts)*exptime
	counts=rebin(xmm_channels,counts,xmmbins)
	
	if errorcalc:
		# Find perturbed spectra for errors
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

newlist=[]
xmm_array=n.array(xmm_datalist)
means=n.mean(xmm_array,axis=0)
meanerrors=[sigma/(2.0*(totalspectra)**0.5) for sigma in n.std(xmm_array,axis=0)]

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
		newrow=[(c-m)/m for c,m in zip(row,means)]
		newlist.append(newrow)

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
			meanval=sum(vals)/len(vals)
			variance=0.
			for val in vals:
				variance+=(val-meanval)**2/float(nerrors-1)
			error=variance**0.5
			errors[i].append(error)

#Transpose the spectrum matrix for processing
xmm_array=n.transpose(n.array(newlist))

print '\nCalculating PCA...'

U,A,V=svd(xmm_array) # This line is actually the entire decomposition
U=n.transpose(U) # U contains all the output components, A contains the variability information
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
	specs.append(spec)

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

if reconstruct:
	# Why are you even reading this far down? This stuff is a mess. Leave while you still can!
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

	if minmax:
		# Find the minimum and maximum spectra for each component
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
				# DO NOT TRUST THESE ERRORS.
				print 'WARNING: All errors on the min/max spectra are extremely dodgy.'
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
		if plot:
			p.show()
		


legendlist=['1st order','2nd order','3rd order','4th order','5th order','6th order','What are you even doing? Nothing has this many orders.','If you were wondering, I was counting zeroth order as being the mean spectrum']

if plot:
	# Seriously, this is a mess as well. It writes the output to some files, just use those and make your own damn plots. They're called comp*_spectrum.txt.
	print '\nPlotting eigenvectors...'

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

	ax2=p.subplot(111)
	
	p.xlabel('Eigenvector',fontsize=axisfontsize)
	p.ylabel('Fractional Variability',fontsize=axisfontsize)
	eigenvalmarker='o'
	for i in range(0,nspec):
		p.plot(i+1,[eigenvals[i]],marker=eigenvalmarker,color=colours[i],ls='-',ms=8)
	p.plot(range(nspec+1,len(eigenvals)),eigenvals[nspec:-1],marker=eigenvalmarker,ls='None',color='y',ms=5)
	ax2.set_yscale('log')
	
	
	formatter=FuncFormatter(to_percent)

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
	# Hmmm. Not touched this for about a year. I wonder what it does.
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


print '\n' # <--- So pointless...





















































































# This is a long way down. I'm scared!










































































# Congratulations! You win ascii Darth!


   #_________________________________
  #|:::::::::::::;;::::::::::::::::::|
  #|:::::::::::'~||~~~``:::::::::::::|
  #|::::::::'   .':     o`:::::::::::|
  #|:::::::' oo | |o  o    ::::::::::|
  #|::::::: 8  .'.'    8 o  :::::::::|
  #|::::::: 8  | |     8    :::::::::|
  #|::::::: _._| |_,...8    :::::::::|
  #|::::::'~--.   .--. `.   `::::::::|
  #|:::::'     =8     ~  \ o ::::::::|
  #|::::'       8._ 88.   \ o::::::::|
  #|:::'   __. ,.ooo~~.    \ o`::::::|
  #|:::   . -. 88`78o/:     \  `:::::|
  #|::'     /. o o \ ::      \88`::::|   "He will join us or die."
  #|:;     o|| 8 8 |d.        `8 `:::|
  #|:.       - ^ ^ -'           `-`::|
  #|::.                          .:::|
  #|:::::.....           ::'     ``::|
  #|::::::::-'`-        88          `|
  #|:::::-'.          -       ::     |
  #|:-~. . .                   :     |
  #| .. .   ..:   o:8      88o       |
  #|. .     :::   8:P     d888. . .  |
  #|.   .   :88   88      888'  . .  |
  #|   o8  d88P . 88   ' d88P   ..   |
  #|  88P  888   d8P   ' 888         |
  #|   8  d88P.'d:8  .- dP~ o8       |   
  #|      888   888    d~ o888    LS |
  #|_________________________________|


#Courtesy of:
#http://www.ascii-art.de/ascii/s/starwars.txt








# Genuine end of the file right here. Awesome going! It's been fun! Have a nice life!