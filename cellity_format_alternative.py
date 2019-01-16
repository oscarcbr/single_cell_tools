#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  cellity_features_format.py
#  
#  Copyright 2019 oscar_bedoya-reina <oscar@oscar-J53kOiSr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


def htseqLFlsToGRM(lFlsHtseq,outGRMFlGns,outGRMFlSpkIns, \
	spkInsPrfx='ERCC-',incldAllGns=True,nmbrSpkIns=2,clcTPM=False, \
	sExcldSmples=None,sIncldSmples=None,rowsNm='file',spltChrc=',', \
	pthCmmnGnNmToENSEMBLG=None,sExcldGnNm=None,excldNmPrfx=None, \
	smplNm=None,nrmlzDpth=False):
	"""
	Input: lFlsHtseq is a list of files for different samples, each of 
	these contains the expression for cells. outGRMFlGns is an output 
	file with all genes expressions in GRM/cellity format. 
	outGRMFlSpkIns is an output file with the expression values for 
	spike-ins in GRM format. spkInsPrfx is the prefix to identify the
	spike-ins. incldAllGns is a switch, if True will include all the 
	genes in the expression results, else is only going to include those
	that are expressed in at least one cell within all samples.
	nmbrSpkIns is the minimal numbero of expressing spike-ins to 
	include the cell. sExcldSmples is a set of samples to exclude from 
	the output. sIncldSmples is a set of samples to include from the 
	output. rowsNm is the name to include with the variable in rows.
	spltChrc is the split character to join the output files.
	pthCmmnGnNmToENSEMBLG is a switch, is True it will transform the
	common gene names to ensembl gene codes. sExcldGnNm is a set of
	gene names to exclude from quantification. excldNmPrfx is a prefix
	to exclude from the cell names. smplNm is the sample name to be
	included if excldNmPrfx is not None, if None is going to include
	the folder parental name. clcTPM is a switch, if True, it will 
	calculate TPMs from RPKMs. nrmlzDpth is a switch, if True, it will
	normalize the counts by read depth.
	Output: outGRMFlGns is an output file with all genes counts in GRM 
	or cellity format. outGRMFlSpkIns is an output file with the counts 
	for spike-ins in GRM format. spkInsPrfx is the prefix to identify 
	the spike-ins.
	NOTE: If excldNmPrfx is not None the sample name is going to be
	added to the header of each sample.
	"""
	if pthCmmnGnNmToENSEMBLG is not None:
		dGnNmENGgnNm = mkDctFrmTb(pthCmmnGnNmToENSEMBLG)
	#
	sSmplCll = set()
	sGnNms = set()
	dSmplCllGnExprsn = {}			
	for htSeqFl in lFlsHtseq:
		if excldNmPrfx is not None:
			if smplNm is None:
				smplNmPrfx = os.path.split(os.path.split(htSeqFl)[-2])[-1]
			else:
				smplNmPrfx = smplNm
			smplCllNm = '.'.join(os.path.split(htSeqFl. \
			replace(excldNmPrfx,''))[1].split('.')[:2])
			smplCllNm = '%s.%s'%(smplNmPrfx,smplCllNm)
		else:
			smplCllNm = '.'.join(os.path.split(htSeqFl)[1].split('.')[:2])
		if sExcldSmples is not None and smplCllNm in sExcldSmples:
			continue
		elif sIncldSmples is not None and smplCllNm not in \
		sIncldSmples:
			continue
		sSmplCll.add(smplCllNm)
		dGnNmExprssn = prseTwoClmnsToDGnNmsVal(htSeqFl,rtrnGrtr0=False, \
		maxKeys=True,dtype=float64)
		if sExcldGnNm is not None:
			dGnNmExprssn = dict([(k,v) for k,v in dGnNmExprssn. \
			items() if k not in sExcldGnNm])
		if pthCmmnGnNmToENSEMBLG is not None:
			dGnNmExprssnTmp,sEnsbmlGnNm = {},set()
			sGnNms = set(dGnNmExprssn.keys())
			while sGnNms:
				k = sGnNms.pop()
				v = dGnNmExprssn.pop(k)
				if dGnNmENGgnNm.has_key(k):
					ensmblGnNm = dGnNmENGgnNm[k]
					if ensmblGnNm in sEnsbmlGnNm:
						print Warning('%s gene has a duplicated name %s'% \
						(k,ensmblGnNm))
					else:
						dGnNmExprssnTmp[ensmblGnNm] = v
				else:
					print Warning \
					('%s gene has no ENSEMBL gene name and was appended as is'% \
					k)
					dGnNmExprssnTmp[k] = v
			dGnNmExprssn = dGnNmExprssnTmp
		dSmplCllGnExprsn[smplCllNm] = dGnNmExprssn
		sGnNms.update(set(dGnNmExprssn.keys()))
	#format input
	ar_srtdSmplCll = array(sorted(sSmplCll))
	ar_srtdGnNms = array(sorted(sGnNms))
	sSpkIns = set([gnNm for gnNm in ar_srtdGnNms if \
	gnNm.find(spkInsPrfx)==0])
	indxSpkIns = array([pos for pos in xrange(len(ar_srtdGnNms)) if \
	ar_srtdGnNms[pos] in sSpkIns])
	#
	ar_gnNmsSmplCll = zeros((len(ar_srtdGnNms),len(ar_srtdSmplCll)))
	lenSrtdGnNms = len(ar_srtdGnNms)
	lenSmplCll = len(ar_srtdSmplCll)
	for gnNmPos in xrange(lenSrtdGnNms):
		gnNm = ar_srtdGnNms[gnNmPos]
		for smplCllPos in xrange(lenSmplCll):
			smplCll = ar_srtdSmplCll[smplCllPos]
			if dSmplCllGnExprsn[smplCll].has_key(gnNm):
				ar_gnNmsSmplCll[gnNmPos][smplCllPos] = \
				dSmplCllGnExprsn[smplCll][gnNm]
	#select only cell with spike-ins expression
	sum_spkInsSmplCll = sum(ma.masked_greater(ar_gnNmsSmplCll[indxSpkIns,:],0).mask,axis=0)
	smplCllsToMsk = ma.masked_greater_equal(sum_spkInsSmplCll,nmbrSpkIns).mask
	#convert to TPM
	if clcTPM:
		ar_sumRPKM = sum(ar_gnNmsSmplCll,axis=0,dtype=float64)
		ar_gnNmsSmplCll = divide(ar_gnNmsSmplCll*1e6,ar_sumRPKM)#TPM
	#normalize by depth
	if nrmlzDpth:
		ar_sumCvrg = sum(ar_gnNmsSmplCll,axis=0,dtype=float64)
		ar_gnNmsSmplCll = divide(ar_gnNmsSmplCll,ar_sumCvrg)
	#write out file
	hdr = '%s%s%s'%(rowsNm,spltChrc,spltChrc.join([ar_srtdSmplCll[p] \
	for p in xrange(len(ar_srtdSmplCll)) if smplCllsToMsk[p]]))
	opndOutGRMFlGns = open(outGRMFlGns,'w')
	opndOutGRMFlSpkIns = open(outGRMFlSpkIns,'w')
	opndOutGRMFlGns.write('%s\n'%hdr)
	for gnNmPos in xrange(lenSrtdGnNms):
		gnNm = ar_srtdGnNms[gnNmPos]
		lExprssns = spltChrc.join([str(ar_gnNmsSmplCll[gnNmPos][p]) for p \
		in xrange(lenSmplCll) if smplCllsToMsk[p]])
		opndOutGRMFlGns.write('%s%s%s\n'%(gnNm,spltChrc,lExprssns))
		if gnNm in sSpkIns:
			opndOutGRMFlSpkIns.write('%s%s%s\n'%(gnNm,spltChrc,lExprssns))	
	#write gene and spike-ins expression 
	opndOutGRMFlGns.close()
	opndOutGRMFlSpkIns.close()
	return 0


def mkDctFrmTb(inTabFl,spltChrc='\t',hdr=True,keyClmn=0,valClmn=1,
	clcStat=None):
	"""
	Method to make dictionaries out of text tables
	Input: inTabFl is a tab file of gene counts with at least two 
	columns. spltChrc is the split character to join the output files. 
	hdr is a switch, if True there is a header and the first row is 
	going to be excluded from the output dictionary. keyClmn is the
	column with the key (0-base). valClmn is the column with the value 
	(0-base). clcStat is a method to calculate an statistic on values
	for each key.
	Output: dKeyVal is dictionary of keys and 1 value.
	NOTE: In case of multiple values and absence of a clcStat method,
	the value return is always going to be the first element after 
	sorting.
	"""
	dKeyVal = {}
	for l in open(inTabFl,'r'):
		if hdr:
			hdr = False
		elif l.strip():
			l = l.splitlines()[0].split(spltChrc)
			k = l[keyClmn]
			v = l[valClmn]
			if dKeyVal.has_key(k):
				dKeyVal[k].append(v)
			else:
				dKeyVal[k]=[v]
	#
	if clcStat is not None:
		dKeyVal = dict([(k,clcStat(v)) for k,v in dKeyVal.items()])
	else:
		dKeyVal = dict([(k,sorted(v)[0]) for k,v in dKeyVal.items()])
	return dKeyVal
	
def prseTwoClmnsToDGnNmsVal(inTsv,gnNmCol=0,valCol=1,rtrnGrtr0=True,
	indxVlCol=False,maxKeys = False,dtype=float32,chrctrSplt=False):
	"""
	Method to parse tsv tables of two columns and make a dictionary of 
	{gnNm:val}.
	Input: inTsv is a two columns table separated by a space/tab. 
	gnNmCol is the column with the key (gene name), valCol is the column
	with the value. rtrnGrtr0 is a switch, if True will only return 
	values greater than 0. indxVlCol is a switch, if str the method will
	select as the valCol number the column in the header with the 
	"indxVlCol" str. avrgKeys is a switch, if True and gene names the 
	maximum expression is going to be taken for genes that are 
	multiplicated. dtype is the type of float to include in parsing the
	data. chrctrSplt is the character to be used to split the lines. 
	Output: dGnNmsVal is a dictionary of keys in gnNmCol and values is 
	valCol
	"""
	#
	dGnNmsVal = {}
	for l in open(inTsv,'r'):
		if indxVlCol and l.split()[0]=='tracking_id':
			valCol = [v for v in l.split() if v.strip()]
			valCol = valCol.index(indxVlCol)
		elif l.strip() and l[0]!='#' and l.split()[0] not in \
		{'tracking_id','gene_id','gene'}:
			if chrctrSplt:
				l=l.splitlines()[0].split(chrctrSplt)
			else:
				l=l.splitlines()[0].split()
			gnNm = l[gnNmCol]
			val = dtype(l[valCol])
			if dGnNmsVal.has_key(gnNm):
				if not maxKeys:
					raise \
					Exception('Key %s is duplicated in the input table' \
					%gnNm)
			elif rtrnGrtr0:
				if val>0:
					if maxKeys:
						if dGnNmsVal.has_key(gnNm):
							if val > dGnNmsVal[gnNm]:
								dGnNmsVal[gnNm] = val
						else:
							dGnNmsVal[gnNm] = val
					else:
						dGnNmsVal[gnNm] = val
			else:
				if maxKeys:
					if dGnNmsVal.has_key(gnNm):
						if val > dGnNmsVal[gnNm]:
							dGnNmsVal[gnNm] = val
					else:
						dGnNmsVal[gnNm] = val
				else:
					dGnNmsVal[gnNm] = val
	return dGnNmsVal
