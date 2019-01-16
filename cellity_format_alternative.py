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

import gzip,os
from numpy import array,float32,float64,mean,sum,std,where,ma,resize, \
zeros,logical_not,divide,inf


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

def clcClltyFeatrs(ar_bamFlPths,ar_fastQflPths,gtfFlPths,clltyFtrsFl, \
	ar_smplNm,ar_cllNm,ar_outFls=None,ovrWrtFls=False):
	"""
	Method to compute the features required for cellity
	Input: ar_bamFlPths is an array of bam files paths for each 
	sample/cell. ar_fastQflPths is an array of fastq paths for each 
	sample/cell. gtfFlPths is a gtf file path. clltyFtrsFl is the output
	file where to write the features results. ar_smplNm is an array of
	sample names. ar_cllNm is an array of cell names. ar_outFls is an
	array of out file names. ovrWrtFls is a switch, if True will 
	overwrite the results files.
	Ouput: clltyFtrsFl is the output file where to write the features 
	results.
	"""
	srtKeys = ['total','mapped','unmapped','unique','multi','intergenic', \
	'intragenic','exonic','intronic','ambigious','exonicM','alignments', \
	'multi-intergenic','multi-intragenic','multi-exonic', \
	'multi-intronic','multi-ambigious','perfect','partly_perfect', \
	'mapped_no_correct','S_0','S_1','S_2','S_3','S_4','S_5','S_6','S_7', \
	'S_8','S_9','S_10+','I','D','INDEL']
	opndClltyFtrsFl = open(clltyFtrsFl,'w')
	opndClltyFtrsFl.write('%s\n'%'\t'.join(['cell','total','mapped', \
	'unmapped','unique','multi','intergenic','intragenic','exonic', \
	'intronic','ambigious','exonicM','alignments','multi.intergenic', \
	'multi.intragenic','multi.exonic','multi.intronic','multi.ambigious', \
	'perfect','partly_perfect','mapped_no_correct','S_0','S_1','S_2', \
	'S_3','S_4','S_5','S_6','S_7','S_8','S_9','S_10.','I','D','INDEL']))
	assert ar_bamFlPths.shape==ar_fastQflPths.shape
	assert ar_smplNm.shape==ar_fastQflPths.shape
	assert ar_cllNm.shape==ar_fastQflPths.shape
	assert ar_outFls.shape==ar_fastQflPths.shape
	nBamFl = ar_bamFlPths.shape[0]
	for posBamFl in xrange(nBamFl):
		bamFl,fastQfl = ar_bamFlPths[posBamFl],ar_fastQflPths[posBamFl]
		smpl,cll = ar_smplNm[posBamFl],ar_cllNm[posBamFl]
		smplCllNm = '%s.%s'%(smpl,cll)
		if ar_outFls is not None:
			outFl = ar_outFls[posBamFl]
			if os.path.exists(outFl) and not ovrWrtFls:
				dStats=mk_dStats(outFl)
			else:
				dStats = generate_mapping_stats(bamFl,gtfFlPths, \
				fastQfl=fastQfl,outFl=outFl,sample_name=str(smpl))
		else:
			dStats = generate_mapping_stats(bamFl,gtfFlPths, \
			fastQfl=fastQfl,outFl=outFl)
		lVls = [str(dStats[k]) for k in srtKeys]
		opndClltyFtrsFl.write('%s\t%s\n'%(smplCllNm,'\t'.join(lVls)))
	opndClltyFtrsFl.close()
	return 0


def generate_mapping_stats(input_file,gtf_file,count=False, \
	geneFeature='transcript',fastQfl=False,sample_name=False, \
	output_file=False,outFl=None,incldERCC=False):
	"""
	Method to generate the mapping stats
	"""
	print("Loading GTF file...") 
	gtf_dict = load_GTF(gtf_file) 
	print("Loaded.")
	#OUTPUT TABLE CONTAING STATS
	output_table = OrderedDict()
	#Dict indicating to which genes a specific read maps to
	#It is a temporary dict
	reads_mapped_to = defaultdict(str)
	exonic_mappings_temp = defaultdict(str)
	#Dict indicating which read is multi-mapped
	#It is a temporary dict
	multi_maps = defaultdict(str)
	sam_files = [input_file]
	exonic_multi_table = defaultdict(str)
	#MAPPABILITY
	output_table["total"] = 0 
	output_table["mapped"] = 0
	output_table["unmapped"] = 0
	output_table["unique"] = 0
	output_table["multi"] = 0
	#CODING VERSUS NON-CODING REGIONS
	output_table["intergenic"] = 0
	output_table["intragenic"] = 0
	output_table["exonic"] = 0
	output_table["intronic"] = 0
	output_table["ambigious"] = 0
	#CODING REGIONS MAPPABILITY
	output_table["exonicU"] = 0
	output_table["exonicM"] = 0
	#ALIGNMENT CODING VS NONCODING
	output_table["alignments"] = 0
	output_table["multi-intergenic"] = 0
	output_table["multi-intragenic"] = 0
	output_table["multi-exonic"]=0
	output_table["multi-intronic"]=0
	output_table["multi-ambigious"]=0
	#ERROR
	output_table["perfect"] = 0
	output_table["partly_perfect"] = 0
	output_table["mapped_no_correct"] = 0
	for i in range(0,10):
		output_table["S_"+str(i)] = 0
	#
	output_table["S_10+"] = 0
	output_table["I"] = 0
	output_table["D"] = 0
	output_table["INDEL"] = 0
	if incldERCC:
		output_table["ercc"] = 0
	#set counters
	reads = Counter()
	multi_reads = defaultdict(str)
	#add to keep track on the aligned reads for total reads values
	if fastQfl:
		sAlgndRds = set()
	#SAM PARSE SAM FILE
	for sam in sam_files:
		print("Parsing bam/sam file...")
		with pysam.AlignmentFile(str(sam),'rb') as f:
			for line in f:
				line = line.tostring()
				split = line.split("\t")
				if (not line.startswith("@PG") and not \
				line.startswith("@HD") and not line.startswith("@SQ") \
				and len(split) >= 10):
					read_name=split[0]
					flagCode = int(split[1])
					chrom = split[2]
					pos = split[3] 
					errors=split[5]
					read=split[9]
					#
					if fastQfl:
						sAlgndRds.add(read_name)
					errors_a = list(errors)
					number = ""
					num = 0
					error_table = defaultdict(int)
					name_and_flag = read_name
					#CHECK IF READ MAPPED OR UNMAPPED 
					#IT US UNMAPPED
					if(flagCode & 0x0004 != 0):
						output_table["unmapped"] += 1
						output_table["total"] += 1
						error_table["*"] += 1
					#IT IS MAPPED
					else:
						if (flagCode & 0x0001 !=0):#This is paired end sequencing
							if (flagCode & 0x0040 != 0):#1st read
								name_and_flag += ";first"
							if (flagCode & 0x0080 != 0):#2nd read
								 name_and_flag  += ";second"
						# CHECK TO WHICH GENE(S) IT MAPPED TO
						genes_info, num_genes, num_exons = get_gene(gtf_dict, [chrom, pos],geneFeature)				   
						#GENE COUNTS: only NON-overlapping genes and if read within exons
						if (count and num_genes == 1 and num_exons > 0):
						  info = genes_info[0]							
						  gene_id = info[4]
						  mapped_to = []
						  if (name_and_flag in reads_mapped_to):
							  mapped_to = reads_mapped_to[name_and_flag]
						  mapped_to.append(gene_id)
						  reads_mapped_to[name_and_flag] = mapped_to
						output_table["alignments"] += 1.0
						#STATS
						if(name_and_flag not in reads):
							reads[name_and_flag] += 1
							output_table["unique"] += 1
							output_table["total"] += 1
							output_table["mapped"] += 1
							if(num_genes == 0):
								output_table["intergenic"] += 1
							elif (num_genes == 1):
								output_table["intragenic"] += 1
								#include ercc reads
								if incldERCC and \
								genes_info[0][4].find('ERCC-')>-1:
									output_table["ercc"] += 1
								#
								if (num_exons==0):
									output_table["intronic"] += 1
								else:
									output_table["exonic"] += 1 
									output_table["exonicU"] += 1
									d = []
									if (name_and_flag in exonic_mappings_temp):
										d = exonic_mappings_temp[name_and_flag]
									d.append([genes_info[0], chrom, pos])
									exonic_mappings_temp[name_and_flag] = d
							elif (num_genes > 1):
								output_table["ambigious"] += 1
						#READ IS MULTI-MAPPED
						else:
							if(reads[name_and_flag] == 1):
								output_table["unique"] -= 1
								output_table["exonicU"] -= 1
								output_table["multi"] += 1
							reads[name_and_flag] += 1
							d = []
							#GET KNOWLEDGE IF FIRST MAPPING EXONIC OR INTRONIC
							if (name_and_flag in exonic_mappings_temp):
								d = exonic_mappings_temp[name_and_flag]
							#output_table["alignments"] += 1.0 
							if(num_genes == 0):
								output_table["multi-intergenic"] +=(1)
							elif (num_genes == 1):
								output_table["multi-intragenic"] += (1)
								if (num_exons==0):
									output_table["multi-intronic"] += (1)
								else:
									output_table["multi-exonic"] += (1) 
									d.append([genes_info[0], chrom, pos])
							elif (num_genes > 1):
								output_table["multi-ambigious"] += (1)
							#IF AT LEAST ONE EXONIC ALIGNMENT
							if (len(d) > 0):
								exonic_multi_table[name_and_flag] = d 
						#PARSE MAPPING ERRORS		 
						for i in errors_a:
							if (re.match("[0-9]",i)):
								number +=(i)
							elif(re.match("[A-Z]",i)):
								num = int(number)
								error_table[i] += num 
								number = ""
						#print output_table
						#TABLE OF HOW MANY READS MAP PERFECT, PARTLY PERFECT, SUBSTITUINTS ETC
						if("M" in  error_table and len(error_table)==1):
							output_table["perfect"] += 1
						elif("M" in error_table and len(error_table) > 1): 
							output_table["partly_perfect"] += 1
						elif("M" not in error_table and "*" not in error_table):
							output_table["mapped_no_correct"] += 1
						if("S" in error_table):
							if(int(error_table["S"]) < 10):
								output_table["S_"+str(error_table["S"])] += 1
							else:
								output_table["S_10+"] += 1
						elif("S" not in error_table):
							output_table["S_0"] += 1

						if("I" in error_table):
							output_table["I"] += 1

						if("D" in error_table):
							output_table["D"] += 1
	
						if("I" in error_table or "D" in error_table):
							 output_table["INDEL"] += 1
	#add total reads	
	if fastQfl:
		sAll = set()
		with pysam.FastxFile(str(fastQfl)) as fh:
			for entry in fh:
				rdNm = entry.name
				sAll.add(rdNm)
				if rdNm not in sAlgndRds:
					output_table["unmapped"] += 1
					output_table["total"] += 1
					error_table["*"] += 1
		#~ print 'len(sAll)->',len(sAll)
	#WEIGHT COUNTS 
	#~ dGnNmCnts = {}
	#~ if (count):
		#~ counts, counts_unique, counts_multi = weight_counts(reads_mapped_to)
		#~ write_counts(output_file+".counts.unique", counts_unique, sample_name)
		#~ write_counts(output_file+".counts", counts, sample_name)
		#~ write_counts(output_file+".counts.multi", counts_multi, sample_name)		
	#~ o = ""
	exonicM = len(exonic_multi_table.keys())
	output_table["exonicM"] = exonicM
	if outFl is not None:
		outFl = str(outFl)
		if not sample_name:
			sample_name = '.'.join(os.path.split(outFl)[1].split('.')[:-1])
		write_stats(outFl,output_table,sample_name)
	dStats = rtrn_dStats(output_table)
	return dStats

def load_GTF(gtf_file):
	"""
	Method for making an interval tree from a gtf file.
	"""
	gtf_index = defaultdict()
	with open(gtf_file) as f:
		for line in f:
			if (not line.startswith("#")):
				entry = line.split("\t")
				entry_addition = entry[8]
				entry_addition = entry_addition.split(";")
				entry_addition = entry_addition[0].split(" ")
				gene_id = entry_addition[1]
				type = entry[2] 
				#TYPE(Gene, exon etc.), START, END, STRAND, gene_ID
				info = [type, entry[3], entry[4], entry[6], gene_id]
				#Build GTF INDEX
				if (type != "" and entry[3]!= entry[4]):
					index = IntervalTree()
					if (entry[0] in gtf_index):
						index = gtf_index[entry[0]]
					index.addi(int(info[1]),int(info[2]),info) 
					gtf_index[entry[0]] = index
	return (gtf_index)

def write_stats(output_file, stats_table, o):
	"""
	Method to write output files
	"""
	#OUTPUT STATS
	f = open(output_file,'w')
	for k,v in stats_table.items():
		if (str(k) in ["unique", "multi", "intragenic", "intergenic", \
		"exonic", "intronic", "ambigious", "exonicM", "exonicU"]):
			v = (v+0.0) / (stats_table["mapped"]+0.0)
			v = '%.2f' % (100.0*(v))
		if (str(k) in ["multi-intragenic", "multi-intergenic", \
		"multi-exonic", "multi-intronic", "multi-ambigious"]):
			v = (v+0.0)
			if (stats_table["alignments"] !=0):
				v = v / (stats_table["alignments"]+0.0)  
			v = '%.2f' % (100.0*(v))
		val = v 
		print str(k) + ":" + str(val)
		o += "," + str(val)
	o += "\n" 
	f.write(o)
	f.close()

def write_counts(output_file, counts,o):
	"""
	Method to write counts to file
	"""
	f = open(output_file,'w')
	o="Gene,"+o
	for k,v in counts.items():
		o = str(k).replace("\"", "")+ "," + str(v) + "\n"
		f.write(o)
	f.close()



def weight_counts(reads_mapped_to):
	"""
	Method weight counts of reads mapped to genes
	"""
	#Considers also multi-mapped reads
	#If multi-mapped, counts are shared across genes
	#Only considers genes that do not have any overlapping genes
	#Considers paired-end as single ended reads
	counts = defaultdict(float)
	counts_unique = defaultdict(float)
	counts_multi = defaultdict(float) 
	#OUTPUT GENE COUNTS
	for read in reads_mapped_to.keys():
		genes = reads_mapped_to[read]
		if (len(genes) == 1):
			counts_unique[genes[0]] += 1.0
			counts[genes[0]] += 1.0
		else:
			for i, g in enumerate(genes):
				counts_multi[g] += 1.0
				counts[g] += 1.0/len(genes)
	return counts, counts_unique, counts_multi


def count_genes(genes_info, reads_mapped_to, name_and_flag):
	"""
	Method to count genes
	"""
	#GENE COUNTS: only NON-overlapping genes
	if (len(genes_info)== 1):
		info = genes_info[0]                            
		gene_id = info[4]
		mapped_to = []
		if (name_and_flag in reads_mapped_to):
			mapped_to = reads_mapped_to[name_and_flag]
		mapped_to.append(gene_id)
		reads_mapped_to[name_and_flag] = mapped_to


def get_gene(gtf_dict,pos_pair,geneFeature='gene'):
	"""
	Method to retrieve genes
	"""
	set_genes = set()
	set_exons = set()
	num_genes = 0
	num_exons = 0
	if (pos_pair[0] not in gtf_dict):
		print ("Ignored pos: "+pos_pair[0])
		return([[], num_genes, num_exons])
	entries = gtf_dict[pos_pair[0]]
	pos = pos_pair[1]
	found = []
	found = entries.search(int(pos_pair[1]))
	list = []
	for e in found:
		info = e[2]
		if(info[0]== geneFeature):
			list.append(info)
			set_genes.add(info[4])
			#~ num_genes+= 1
		elif (info[0] == "exon"):
			#~ num_exons+=1
			set_exons.add((info[1],info[2]))
	num_genes = len(set_genes)
	num_exons = len(set_exons)
	return([list, num_genes, num_exons])


def rtrn_dStats(stats_table):
	"""
	Method to return a dictionary of stats
	"""
	dStats = {}
	for k,v in stats_table.items():
		if (str(k) in ["unique", "multi", "intragenic", "intergenic", "exonic", "intronic", "ambigious", "exonicM", "exonicU"]):
			v = (v+0.0) / (stats_table["mapped"]+0.0)
			v = '%.2f' % (100.0*(v))
		if (str(k) in ["multi-intragenic", "multi-intergenic", "multi-exonic", "multi-intronic", "multi-ambigious"]):
			v = (v+0.0)
			if (stats_table["alignments"] !=0):
				v = v / (stats_table["alignments"]+0.0)  
			v = '%.2f' % (100.0*(v))
		val = v 
		dStats[k]=v
	return dStats
