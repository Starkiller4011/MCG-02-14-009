#!/bin/bash

tstep=$1
outdir=$2
# table=$3
obsid=$3
source=$4

table=${obsid}/EPICclean.fits

# source=ngc1365
# 
export SAS_CCF=/scratch/mlparker/${source}/${obsid}/ccf.cif
# export SAS_ODF=/scratch/mlparker/${source}/0029740101/0301_0029740101_SCX00000SUM.SAS

# xmm_unfiltered_gti=/scratch/mlparker/${source}/0029740101/0301_0029740101_EPN_S003_ImagingEvts.ds

# table=/scratch/mlparker/${source}/0029740101/pn/0301_0029740101_EPN_S003_ImagingEvts.ds
# table=/scratch/mlparker/${source}/0029740701/pn/0302_0029740701_EPN_S003_ImagingEvts.ds
# table=/scratch/mlparker/${source}/0029740801/pn/0303_0029740801_EPN_S003_ImagingEvts.ds

count=0

# tstep=20000

# outdir=products_pca_8

let nfiles=$(ls -1 ${outdir}/${source}_pn_tstep${tstep}_no*.gti | wc -l)-1
for k in $(seq 0 $nfiles); do
	outgti=${outdir}/${source}_pn_tstep${tstep}_no${count}.gti
	
	srcreg=`tail -1 ${obsid}/src.reg | head -1`
	bkgreg=`tail -1 ${obsid}/bkg.reg | head -1`
	
	echo -e '\nExtracting source spectrum...'
	# Extract source spectrum
	evselect table=$table \
		withspectrumset=yes \
		spectrumset=${outdir}/${source}_pn_tstep${tstep}_no${count}_src.fits \
		energycolumn=PI \
		spectralbinsize=5 \
		withspecranges=yes \
		specchannelmin=0 \
		specchannelmax=20479 \
		expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN '${srcreg}') && gti('$outgti',TIME)'
	
	echo -e '\nExtracting background...'
	#Extract background spectrum
	evselect table=$table \
		withspectrumset=yes \
		spectrumset=${outdir}/${source}_pn_tstep${tstep}_no${count}_bkg.fits \
		energycolumn=PI \
		spectralbinsize=5 \
		withspecranges=yes \
		specchannelmin=0 \
		specchannelmax=20479 \
		expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN '${bkgreg}') && gti('$outgti',TIME)'
	
	backscale spectrumset=${outdir}/${source}_pn_tstep${tstep}_no${count}_src.fits \
		badpixlocation=$table
	
	backscale spectrumset=${outdir}/${source}_pn_tstep${tstep}_no${count}_bkg.fits \
		badpixlocation=$table
		
	echo -e '\nGenerating response files...'
	rmfgen spectrumset=${outdir}/${source}_pn_tstep${tstep}_no${count}_src.fits \
		rmfset=${outdir}/${source}_pn_tstep${tstep}_no${count}.rmf
	
	arfgen spectrumset=${outdir}/${source}_pn_tstep${tstep}_no${count}_src.fits \
		arfset=${outdir}/${source}_pn_tstep${tstep}_no${count}.arf \
		withrmfset=yes \
		rmfset=${outdir}/${source}_pn_tstep${tstep}_no${count}.rmf \
		badpixlocation=$table \
		detmaptype=psf
	
	echo -e '\nGrouping files with pha...'
	grppha ${outdir}/${source}_pn_tstep${tstep}_no${count}_src.fits ${outdir}/${source}_pn_tstep${tstep}_no${count}_src.pha <<EOF
chkey BACKFILE ${source}_pn_tstep${tstep}_no${count}_bkg.fits
chkey RESPFILE ${source}_pn_tstep${tstep}_no${count}.rmf
chkey ANCRFILE ${source}_pn_tstep${tstep}_no${count}.arf
write !${outdir}/${source}_pn_tstep${tstep}_no${count}_src.pha
exit !${outdir}/${source}_pn_tstep${tstep}_no${count}_src.pha
EOF

# 
# 	rm ${outdir}/${source}_pn_tstep${tstep}_no${count}.fits
		#rm xmm/products_tstep/PN_spec_${count}.fits
		#mv xmm/products_tstep/PN_spec_${count}_out.fits products_tstep/PN_spec_${count}.fits
		
	let count=$count+1
		
done
