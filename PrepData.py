#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import RFE
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import normalize
from sklearn.preprocessing import Imputer
import zlib,bz2
import sys
from math import log
import itertools
import seaborn as sns
from itertools import combinations
from itertools import product
from itertools import combinations_with_replacement
from itertools import islice
from collections import OrderedDict
from Bio.Seq import Seq
from itertools import islice
import statsmodels.api as sm
from statsmodels.formula.api import ols
from collections import defaultdict
import MySQLdb as mdb
import glob
import re

def reverseSeq(seq):
        orig=Seq(seq)
        rev=orig.reverse_complement()
        return rev

def window(seq, n):
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result =result[1:] + (elem,)
        yield result



def trifonov(string):
	prod=1e9	##distribution looks normal/gaussian, but this scale gives better 'looking' numbers. 
			##e.g. from 1.248e-9 to 1.24. All numbers end up as e-9, rescaled to >0 and < 10.
	for n in range(2,len(string)+1):
		if n >=18: argmin=4
		else:	argmin=min(n,3)
			
		cd=set(map(lambda q: "".join(q),list(window(string,n))))
		gg=list(map(lambda q: "".join(q),list(window(string,argmin))))
		
		prod*=len(cd)/float(len(gg))
	return prod

def entropy(sequence):
	lsq=list(window(sequence,2))
	dNs=["".join(item) for item in lsq]
	uniquedN=list(set(dNs))
	if len(uniquedN) <= 1:
		return 0
	probs = [dNs.count(item)/float(len(lsq)) for item in uniquedN] 
	ent = 0.
	# Compute standard entropy.
	for i in probs:
		ent -= i * log(i,len(lsq))
	return ent

def compressString(string):
	clen=len(zlib.compress(string))
	return clen



cnx=mdb.connect('cmmrdbdev.research.sickkids.ca','greg','t43tt43t','gRNA_table')
cur=cnx.cursor()

def Gidfind(Gid):
	try:
		cur.execute("""Select * from gRNA_table.GRNA_PARTITIONS where AI_GRNA = %s""" % (Gid,))
		bits=cur.fetchall()
		bits[0]
	except (mdb.Error,IndexError),e:
		return
	return bits[0]

def GeneFind(marker):
	try:
		cur.execute("""Select gMGI_ID,gEnsemblID,gChr,gGenomicSite from musmusculims.genes where gMGI_Symbol = %s""" % ("\""+marker+"\"",))
		bits=cur.fetchall()
		bits[0]
	except (mdb.Error,IndexError),e:
		return
	return bits[0]
	


def get_crisps(sequence):
	if sequence.strip().startswith("CC"):
		#result.append("R")
		seq=sequence.strip()[3:]
	elif re.match(r'[ACGTN]{21}[G]{2}',sequence.strip()):
		seq=sequence.strip()[:20]
	else:
		return []
		
	finder="./Data/*"+seq+".txt"
	files=glob.glob(finder)
	p=[]
	for f in files:
		io=open(f,'r').readlines()
		try:
			io.pop(0)
			p+=map(lambda l: l.strip(),io)
		except IndexError:
			pass
	found=list(set(p))
	return found



casBase="/var/www/CRISPR-Analyser-master/bin/"
crispFiles=glob.glob("./Data/*.txt")
crispFiles.sort()

data=pd.read_csv("cripsr_data_04_04_17.csv")
def InitClean(df):
	df=df[(df.go_screened > 0)]								
	df=df[(df.secondary_allele == "Deletion") | (df.primary_allele == "Deletion")]	##looking for deletions only
	df=df[(df.experimental) != "t"]		#no experimental data
	print len(df)
	generalColumns=["production_centre",
			"marker_symbol",
			"primary_allele",
			"secondary_allele",
			"delivery_method",
			"cas9_d10a",
			"mrna_protein",
			"mrna_nuclease_conc",
			"grna_conc",
			"grnas_per_gene",
			"grna_sequences",
			"grna_coordinates",
			"grna_cleave_sites",
			"mi_date",
			"g0_obtained_date",
			"genotype_confirmed_date",
			"status"
			]
	countColumns=["e_injected",
			"embryos_survived",	
			"e_transferred",
			"go_screened",
			"total_count_of_mutagenized_g0_detected_by_screen",
			"g0_with_nhej_mutation",
			"g0_with_deletion_mutation",
			"g0_bred",
			"genotype_confirmed_f1s"	
			]			

	df=df[generalColumns+countColumns]	##include both
	for col in countColumns:	df[col].fillna(0,inplace=True)
	rawdata=df.columns.values.tolist()
	df.to_csv("./crispr.clean.csv",sep=";")##We already have comma-separated columns
	#return df




############
#
#######################
#
#
########################
#
############


def checkValid(grna_seqs,Midcoord,chrom):
	valid=[]
	lcl=1
	for Gid in grna_seqs:
			##(84437773L, 5L, 'chr5_0_137957977', 'CCTGCCGCACATCCTGTCACCAA')
		Ginfo=Gidfind(Gid)
		if Ginfo:
			Lid,dummy,location,sequence=Ginfo
			chromosome,strand,cutsite=location.split("_")
			chromosome=chromosome.lstrip("chr")
			if chromosome == "98":
				chromosome="X"
			elif chromosome == "99":
				chromosome = "Y"

			if chrom == chromosome and abs(int(cutsite)-Midcoord) < 400000:
				valid.append([sequence,cutsite])
	
		lcl+=1
	return valid


def grnaID(filename):
	Ic=open(filename,'r').readlines()
	header=Ic[0].strip()
	header+=";chromosome\n"


	fheader=header.split(";")
	fheader[-1]=fheader[-1].strip()

	Ic.pop(0)
	kl=open('crispr.clean.verified.csv','w')
	kl.write(header)

	iCols=["marker_symbol",	
		"grna_sequences",
		"grna_coordinates"
		]
	IndexCol=dict(zip(iCols,[fheader.index(k) for k in iCols]))

	for r in Ic:
		data=r.split(";")
		grna_seqs=map(lambda p: p.strip(),data[IndexCol['grna_sequences']].split(","))
		marker_symbol=data[IndexCol['marker_symbol']].strip()
		grna_seqs.sort()
		grna_ids=[]
		gene_info=GeneFind(marker_symbol)
		if not gene_info:
			pass
		else:
			gMGI_ID,gEnsemblID,gChr,gGenomicSite=gene_info
			seq_positions=map(lambda p: int(p),gGenomicSite.split("-"))
			seq_positions.sort()
			start,end=seq_positions[0],seq_positions[1]
			midpoint=np.mean([start,end])
			
			for s in grna_seqs:
				grna_ids+=get_crisps(s)	
			valids=checkValid(grna_ids,midpoint,gChr)	
			finalseq=[]
			finalcut=[]
			for v in valids:
				mS,mC=v
				vS=reverseSeq(mS)
				if mS in grna_seqs:
					finalseq.append(mS)
					finalcut.append(mC)
				elif vS in grna_seqs:
					finalseq.append(mS)
					finalcut.append(mC)
			if len(grna_seqs) == len(set(finalseq)):
				line_only=r.strip()
				kl.write(line_only+";"+str(gChr)+"\n")
	kl.close()





def cutSizeDet(sorted_cuts):
	cutsize=[]
	if len(sorted_cuts) == 2:
		cutsize=abs(sorted_cuts[0]-sorted_cuts[1])
		cuttype="2"
	elif len(sorted_cuts) == 4:
		Lft=[sorted_cuts[0],sorted_cuts[1]]
		Rgt=[sorted_cuts[2],sorted_cuts[3]]
		cutOuter=abs(sorted_cuts[3]-sorted_cuts[0])
		cutInner=abs(sorted_cuts[2]-sorted_cuts[1])
		###
		cutsize=np.mean([cutOuter,cutInner])
		cuttype="4"
	elif len(sorted_cuts) > 2:
		#print len(sorted_cuts),sorted_cuts
		leftD=[]
		leftC=[]
		rightC=[]
		rightD=[]
		for j in range(0,(len(sorted_cuts)/2)+1):
			leftD.append(sorted_cuts[j+1]-sorted_cuts[j])
			leftC.append(str(sorted_cuts[j])+"_"+str(sorted_cuts[j+1]))
		for j in reversed(range(len(sorted_cuts)/2+1,len(sorted_cuts))):
			rightD.append(sorted_cuts[j]-sorted_cuts[j-1])
			rightC.append(str(sorted_cuts[j-1])+"_"+str(sorted_cuts[j]))

		jntD=leftD+rightD
		jntC=leftC+rightC
		mDist=max(jntD)
		
		for a,b in zip(jntD,jntC):
			if a == mDist:
				flankLeft,flankRight=map(lambda i: int(i),b.split("_"))
				break
		Li=sorted_cuts.index(flankLeft)
		Ri=sorted_cuts.index(flankRight)
		Lft=sorted_cuts[:(Li+1)]
		Rgt=sorted_cuts[Ri:]
		cutOuter=abs(Rgt[-1]-Lft[0])
		cutInner=abs(flankRight - flankLeft)
		###
		cutsize=np.mean([cutOuter,cutInner])
		cuttype="other"
	else:
		print 20*"\tWTF\n"

	return cutsize,cuttype
		
def AggregateExp(filename):
	io=open(filename,'r').readlines()
	header=io[0].split(';')
	header[-1]=header[-1].strip()	#strip from lambda will remove first header which is null,messing up index
	header.append('DeletionSize')
	header.append('DeletionType')
	header.append('AggregExp')
	io.pop(0)
	experiments={}
	Uexperiments={}
	aggregateColumns=["e_injected",
			"embryos_survived",	
			"e_transferred",
			"go_screened",
			"total_count_of_mutagenized_g0_detected_by_screen",
			"g0_with_nhej_mutation",
			"g0_with_deletion_mutation",
			"g0_bred",
			"genotype_confirmed_f1s"	
			]	
	aggIndex=dict(zip(aggregateColumns,[header.index(k) for k in aggregateColumns]))
	Indexagg=dict(zip([header.index(k) for k in aggregateColumns],aggregateColumns))

	iCols=["marker_symbol",	
		"grna_sequences",
		"grna_coordinates",
		"production_centre",
		"grna_cleave_sites"
		]
	IndexCol=dict(zip(iCols,[header.index(k) for k in iCols]))


	aggData={}
	for line in io:
		d=line.strip().split(";")
		grna_seqs=map(lambda p: p.strip(),d[IndexCol['grna_sequences']].split(","))
		cut_sites=map(lambda p: int(p),d[IndexCol['grna_cleave_sites']].split(","))
		cut_sites.sort()
		###create a unique identifier based on Gene,gRNA_cut_sites,production_centre
		gene=d[IndexCol['marker_symbol']].upper()
		uniq="_".join(map(lambda s: str(s),cut_sites))	
		prod_centre="_".join(d[1].split())
		#experiments[d[2]]=uniq
		uniqueID=gene+"_"+prod_centre+"_"+uniq
		cutS,cutT=cutSizeDet(cut_sites)	
		_line=[]
		if uniqueID not in aggData:
			for i,datapoint in enumerate(d):
				if i in Indexagg:	###Is this numeric value??
					_line.append(float(datapoint))
				else:	###
					_line.append(datapoint)
			if cutS > 100000:
				cutS=7000
				cutT="Capped"
			_line.append(cutS)
			_line.append(cutT)
			_line.append(1)
			aggData[uniqueID]=_line
		else:
			for i,datapoint in enumerate(d):
				if i in Indexagg:
					aggData[uniqueID][i]+=float(datapoint)
				else:
					pass	##BELOW WE CAN CHECK IF THERE ARE OTHER EXPERIMENTAL ALTERNATES
					#if i != 0 and datapoint != aggData[uniqueID][i]:
					#	print "\t",gene,i,header[i],datapoint,"\t\t",aggData[uniqueID][i]
				aggData[uniqueID][-1]+=1
	maindat=[]
#	maindat.append(header)
	for k,o in aggData.iteritems():
		maindat.append(o)
	header[0]="Index"
	df=pd.DataFrame(maindat,columns=header)
	df=df[
		(df.production_centre == "JAX") | 
		(df.production_centre == "UCD") |
		(df.production_centre == "WTSI") |
		(df.production_centre == "BCM") |
		(df.production_centre == "TCP") |
		(df.production_centre == "MARC") |
		(df.production_centre == "Harwell") |
		(df.production_centre == "RIKEN BRC") |
		(df.production_centre == "ICS") 
		]
	
	df.to_pickle('crispr.clean.verified.aggregated.DataFrame.pkl')



if __name__ == "__main__":
	InitClean(data)
	grnaID('crispr.clean.csv')
	AggregateExp('crispr.clean.verified.csv')
