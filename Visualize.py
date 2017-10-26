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

#circle=Circle((horizC,lclMean),0.05,color='k')
#patches.append(circle)



data=pd.read_pickle("Pared.DataFrame.nan.pkl")
FounderRate=data.g0_with_deletion_mutation/data.go_screened
data["FounderRate"]=FounderRate
#data=data[(data.DeletionSize < 5000)]
Sbins=[0,50,100,150,200,250,500,1000,2000,3000,5000]
Fbins=list(np.arange(0,.9,0.1))


bin_edges=map(lambda i: float(i),stats.mstats.mquantiles(data.DeletionSize,[np.arange(0,1,0.1)])[0])
bin_edges[0]=0

f_edges=map(lambda i: float(i),stats.mstats.mquantiles(data.FounderRate,[np.arange(0,1,0.1)])[0])
f_edges[0]=0
#print f_edges
#f_edges=list(np.arange(0,1.1,0.1))

def plotOne(idata):
	#ax=sns.regplot(y="FounderRate",x="DeletionSize",data=data,x_estimator=np.mean,truncate=True)
	#plt.show()
	#ax=sns.jointplot("DeletionSize","FounderRate",data=data,kind='kde',space=0,color='b')
	#ax=(sns.jointplot("DeletionSize","FounderRate",data=data,space=0,color='b')).plot_joint(sns.kdeplot,zorder=0,n_levels=8)
	sns.set(style="whitegrid", palette="muted")
	#'delivery_method'
	#'DeletionType'
	#sns.swarmplot(x="production_centre",y="FounderRate",hue="DeletionType",data=data)
	ax=sns.swarmplot(x="production_centre",y="DeletionSize",hue="DeletionType",data=data)
	axes = ax.axes
	axes.set_ylim(0,5000)
	#axes[0,0].set_ylim(0,)
	#axes[0,1].set_ylim(0,)
	plt.show()



def whichBin(value,bins,Btype):
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
				return 0
			elif k+1 == len(bins) or value > cVal:
				return mBin
			elif value > bins[k] and value <= bins[k+1]:
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

def BinItems(df,_S,_F):
	SizeBins={0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}
	ScoreBins={0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}

	for j,k in zip(df.DeletionSize,df.FounderRate):
		Sbin=whichBin(j,_S,"Score")
		Fbin=whichBin(k,_F,"Del")
		SizeBins[Sbin].append(k)
		ScoreBins[Fbin].append(j)
			#ScoreBins[Lbin].append(a)
	return SizeBins,ScoreBins
#sb,Db=BinItems(data,bin_edges,f_edges)
sb,Db=BinItems(data,Sbins,Fbins)

def BinFounders(df,_S):
	SizeBins={0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}
	ScoreBins={0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}
	print len(df)
	df=df[(df.status == "Founders Obtained") | (df.status == "Genotype confirmed")]
	print len(df)
	for j,k in zip(df.DeletionSize,df.FounderRate):
		Sbin=whichBin(j,_S,"Score")
		SizeBins[Sbin].append(k)
	return SizeBins
Qb=BinFounders(data,Sbins)
#print sb
#print Qb
#sys.exit()

#print "\n\n"
counts=[]
edges=[]
FounderR=[]
ConCount=[]
for k,j in sb.iteritems():
	FounderR.append(np.mean(j))
	counts.append(len(j))
	ConCount.append(len(Qb[k]))
	edges.append(Sbins[k])


#print counts
#print ConCount
#sys.exit()
#########################################
###########################################
ScoresR=[]
for k,j in Db.iteritems():
	if not len(j):
		break
	ScoresR.append(np.mean(j))

#FounderR=map(lambda k: k*max(ScoresR),FounderR)

#print ScoresR
#print FounderR


#	print k,len(j),Fbins[k],np.mean(j)
#############################################
##########################################
blues=['#dae8f5', '#bad6ea', '#88bedc', '#539dcc', '#2a7ab9', '#0b559f']	# generated from sns.color_palette("Blues")
reds=['#fddbcb', '#fcaf93', '#fb8161', '#f44e38', '#d52221', '#a91016']		# generated from sns.color_palette("Reds")
blues.reverse()


#fig6 = plt.figure(figsize=(60,800))
fig6 = plt.figure(figsize=(60,80))
#fig6 = plt.figure(figsize=(80,80))
#fig6 = plt.figure()
ax6 = fig6.add_subplot(111)

allvals=[]			
rectangle=[]


#####Create data frame with new values


#
#  fig, ax = plt.subplots()
#    sb.regplot(x='round', y='money', data=firm, ax=ax)
#    ax2 = ax.twinx()
#    sb.regplot(x='round', y='dead', data=firm, ax=ax2, color='r')
#    sb.plt.show()


#normV=map(lambda p: p*max(counts)/10, FounderR)


lclCol=blues[2]
lclF=reds[1]
verticalCount=1
keepDates=[]
horizCount=0
bottomV=-5
circles=[]

for i,c in enumerate(counts):
	rectangle.append(patches.Rectangle((horizCount,0),30,c,facecolor=lclCol,ec='w',linewidth=0.5))

	rectangle.append(patches.Rectangle((horizCount+1,0),28,ConCount[i],facecolor=lclF,ec='w',linewidth=0.5))
	### Annotation ON TOP
	fontsize=max(12,c/15)
	half_inc=15
	#plt.text(horizCount+half_inc,c+max(5,c*.05),str(c),color='k',fontsize=fontsize,ha='center',va='center')
	plt.text(horizCount+half_inc,c+7,str(c),color='k',fontsize=fontsize,ha='center',va='center')
	plt.text(horizCount+half_inc,7,str(ConCount[i]),color='w',fontsize=fontsize,ha='center',va='center')
	#######


	plt.text(horizCount,bottomV,str(edges[i]),fontsize=15,rotation=90,ha='center',va='top')
	#plt.text(horizCount,bottomV,str(edges[i]),fontsize=15,rotation=90,ha='center',va='top')



	plt.plot([horizCount+.05,horizCount+.95],[bottomV,bottomV],linewidth=4,color='k')
	plt.plot([horizCount,horizCount],[0,bottomV],linewidth=3,color='k')

	##
	#circles.append(patches.Circle((horizCount+.5,normV[i]),20,color='g',ec='k',linewidth=.5,zorder=100,transform=None))
	#circles.append(patches.Circle((horizCount+.5,FounderR[i]),2,color='g',ec='k',linewidth=1,zorder=1000))
	#try:
	#	circles.append(patches.Circle((horizCount+.5,ScoresR[i]),2,color='r',ec='k',linewidth=1,zorder=1000))
	#except:
	#	pass

	horizCount+=30
	rowB=[]

#plt.text(horizCount/float(2),-35,"Intended Deletion Size (bps)",fontsize=20,ha='center',va='center')
#plt.show()
#sys.exit()
#for c in circles:
#	ax6.add_patch(c)
for r in rectangle:
        ax6.add_patch(r)
plt.axis('equal')
plt.axis('off')
ax6.set_aspect('auto')
#plt.axis('tight')
plt.show()
#fig6.show()
#fig6.savefig('rect6.png', dpi=90, bbox_inches='tight')
