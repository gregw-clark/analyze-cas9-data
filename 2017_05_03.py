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
import operator
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as patches


circle=Circle((horizC,lclMean),0.05,color='k')
patches.append(circle)

via=open('viabilityReport.csv','r').readlines()
via.pop(0)

gdict={}
for d in via:
	#'Rhox13', 'MGI:1920864', 'MUAA', 'both', 'homozygote', 'Viable'
	data=map(lambda l: l.strip("\""),d.strip().split(","))
	if len(data) == 7:
		marker,mgi_id,dummy,sex,genotype,pheno,conflict=data
	elif len(data) == 6:
		marker,mgi_id,dummy,sex,genotype,pheno=data
	gdict[marker.strip().upper()]=pheno

def getVia(gdict,marker):
	uMark=marker.strip().upper()
	try:
		gdict[uMark]
		return gdict[uMark]
	except KeyError:
		return None


data=pd.read_pickle("Pared.DataFrame.pkl")
#print len(data)
#print data.genotype_confirmed_f1s
#data=data[(data.genotype_confirmed_f1s >= 0) | (data.g0_bred >= 0)]
#print len(data)
#print data.g0_bred
bin_edges = stats.mstats.mquantiles(data["DeletionSize"], [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1,])
aEdges=[int(b) for b in bin_edges]
Phenotype=map(lambda l: getVia(gdict,l),data.marker_symbol)
data["Phenotype"]=Phenotype
#Phenotype=map(lambda l: getVia(gdict,l),data.marker_symbol)
#print Phenotype

GLTRate=data.genotype_confirmed_f1s/data.g0_bred
data["GLTRate"]=GLTRate
FounderRate=data.g0_with_deletion_mutation/data.go_screened
data["FounderRate"]=GLTRate


def binbin(bins,value):
#	print len(bins)
#	sys.exit()
	for k in range(len(bins)):
		if k + 1 == len(bins):
			return len(bins) -1
		elif value >= bins[k] and value < bins[k+1]:
			return k
uD=defaultdict(list)
for z,k in zip(list(GLTRate),list(data["DeletionSize"])):
	if np.isnan(z) == False and np.isinf(z) == False:
		abin=binbin(bin_edges,k)
		uD[abin].append(z)
km=uD.keys()
km.sort()
lk=[]
lsu=[]
lsl=[]
lsm=[]
sizes=[]
for k in km:
		print round(k,3),"\t\t",len(uD[k]),"\t\t",round(np.mean(uD[k]),3)#,len(uD[k])
		lk.append(k)
		ups=np.percentile(uD[k],90)
		downs=np.percentile(uD[k],10)
		lsu.append(ups)
		lsl.append(downs)
		lsm.append(np.mean(uD[k]))
		sizes.append(len(uD[k]))
#for k in range(len(lk)):
#	plt.plot([lk[k],lk[k]],[lsu[k],lsl[k]],linewidth=2,color='k')
#plt.scatter(lk,lsm,s=sizes,color='k')
#plt.show()
#print len(list(GLTRate))
#print len(list(data["DeletionSize"]))
#sys.exit()

Sbins=[0,50,100,150,200,250,500,1000,2000,3000,5000]
Fbins=list(np.arange(0,.9,0.1))

newticks={}
for g in range(len(Sbins)):
	if g == len(Sbins)-1:
		newticks[g]="[5000+)"
	else:
		newticks[g]="["+str(Sbins[g])+","+str(Sbins[g+1])+")"


def whichBin(value,bins,Btype,newticks):
	if Btype=="Score":
		cVal=15000
		mBin=10
	elif Btype=="Del":
		cVal=0.9
		mBin=8
	if Btype == "Score":
		cVal=15000	
		for k in range(len(bins)):
			if value == 0:
				#return newticks[0]
				return 0
			elif k+1 == len(bins) or value > cVal:
				#return newticks[mBin]
				return mBin
			elif value > bins[k] and value <= bins[k+1]:
				#return newticks[k]
				return k
	if Btype == "Del":
		cVal=0.8	
		mBin=8
		for k in range(len(bins)):
			if value == 0:
				return 0
			elif value == 1:
				return mBin
			elif k+1 == len(bins) or value > cVal:
				return mBin
			elif value > bins[k] and value <= bins[k+1]:
				return k

#FounderRateBin=map(lambda k: whichBin(k,Fbins,"Del",newticks),data.FounderRate)
GLTRateBin=map(lambda k: whichBin(k,Fbins,"Del",newticks),data.GLTRate)
DelSizeBin=map(lambda k: whichBin(k,Sbins,"Score",newticks),data.DeletionSize)

#data["FounderRateBin"]=FounderRateBin
data["GLTRateBin"]=GLTRateBin
data["DelSizeBin"]=DelSizeBin



#Sbins=[0,50,100,150,200,250,500,1000,2000,3000,5000]
newticks=[]
for g in range(len(Sbins)):
	if g == len(Sbins)-1:
		newticks.append("[5000+]")
	else:
		newticks.append("["+str(Sbins[g])+","+str(Sbins[g+1])+")")
altticks=[]
for g in range(len(Fbins)):
	if g == len(Fbins)-1:
		altticks.append("[0.8+]")
	else:
		altticks.append("["+str(Fbins[g])+","+str(Fbins[g+1])+")")

#data=data[(data.delivery_method=="Cytoplasmic Injection") | (data.delivery_method=="Pronuclear Injection")]
#ax=sns.boxplot(x="production_centre",y="FounderRate",hue="delivery_method",data=data,palette="muted",color='g')
#ax=sns.swarmplot(x="DelSizeBin",y="FounderRate",data=data)
#ax.set(xticklabels=newticks)
#plt.scatter(data["GLTRate"],data["DeletionSize"])
#plt.show()
#ax=sns.jointplot(x="GLTRate",y="DeletionSize",data=data,kind="kde")
#plt.show()
#sys.exit()

#ax=sns.boxplot(x="DelSizeBin",y="GLTRate",data=data,palette="muted")
#ax=sns.swarmplot(x="DelSizeBin",y="GLTRate",data=data)
#ax.set(xticklabels=altticks)

#ax=sns.boxplot(x="GLTRateBin",y="DeletionSize",data=data,palette="muted")
#ax=sns.swarmplot(x="GLTRateBin",y="DeletionSize",data=data)
#ax.set(xticklabels=altticks)

#plt.xticks(newticks)
ax=sns.boxplot(x="Phenotype",y="FounderRate",data=data,palette="muted")
ax=sns.swarmplot(x="Phenotype",y="FounderRate",data=data)
ax.set_xticks(np.arange(1,4,1))
plt.show()



#for e in data.itertuples():
#	print e.marker_symbol#,e.DeletionSize,e.DelSizeBin
#print DelSizeBin


