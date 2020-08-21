#!/bin/python

import pyfits
import subprocess
import os

timebinsize=10000

source='mrk1048'

#obsids=['0150470601','0690870101','0690870501']
obsids=['0150470601']

directory='/scratch/mlparker/%s/' % source
tabname='/EPICclean.fits'

tables=[]
gtis=[]
for obsid in obsids:
	tables.append(directory+obsid+tabname)
	gtis.append(directory+obsid+tabname)


#tables=['/scratch/mlparker/1H0707/0653510301/EPICclean.fits',\
	#'/scratch/mlparker/1H0707/0653510401/EPICclean.fits',\
	#'/scratch/mlparker/1H0707/0653510501/EPICclean.fits',\
	#'/scratch/mlparker/1H0707/0653510601/EPICclean.fits']

for j in range(0,len(obsids)):
	obsid=obsids[j]
	
	outdir='products_pca_%s/' % obsid
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	outdir+='%ss/' % timebinsize
	if not os.path.exists(outdir):
		os.mkdir(outdir)

	table=tables[j]
	xmm_gti=pyfits.open(gtis[j])
	
	prihdu=xmm_gti[0]
	try:
		gtidata=xmm_gti['STDGTI04'].data
	except:
		gtidata=xmm_gti['STDGTI'].data
	starttime=gtidata[0][0]
	endtime=gtidata[-1][1]

	gtiheader=xmm_gti[6].header

	newstarts=[]
	newends=[]


	#for row in gtidata:
	nbins=int((gtidata[-1][1]-gtidata[0][0])/timebinsize)
	if (gtidata[0][1]-gtidata[-1][0])/timebinsize>float(nbins):
		nbins+=1
	#if (row[1]-row[0])/timebinsize>nbins:
		#nbins+=1
	for i in range(0,nbins):
		newstarts.append(gtidata[0][0]+timebinsize*(i))
		newends.append(gtidata[0][0]+timebinsize*(i+1))

	#print newstarts
	#print newends

	for i in range(0,len(newstarts)):
		outfname=outdir+'%s_pn_tstep%s_no%s.gti' % (source,timebinsize,i)
		if not os.path.exists(outfname):
			pyfits.writeto(outfname,prihdu.data,prihdu.header)
			start=[newstarts[i]]
			stop=[newends[i]]
			
			col1 = pyfits.Column(name='START',format='1D',array=start)
			col2 = pyfits.Column(name='STOP',format='1D',array=stop)
			cols = pyfits.ColDefs([col1,col2])
			tbhdu = pyfits.new_table(cols)
			
			pyfits.append(outfname,tbhdu.data,gtiheader)
	subprocess.call(['./xmm_pca.sh',str(timebinsize),outdir[:-1],obsid,source])
