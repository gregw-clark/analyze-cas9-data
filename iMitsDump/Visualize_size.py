#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from math import log
import itertools
from collections import OrderedDict
from collections import defaultdict
import operator
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from altfuncs import BinItems,BinFounders 
import seaborn as sns

def PlotSizeEfficiencies():
	blues=['#dae8f5', '#bad6ea', '#88bedc', '#539dcc', '#2a7ab9', '#0b559f']	# generated from sns.color_palette("Blues")
	reds=['#fddbcb', '#fcaf93', '#fb8161', '#f44e38', '#d52221', '#a91016']		# generated from sns.color_palette("Reds")
	blues.reverse()

	fig6 = plt.figure(figsize=(60,80))
	ax6 = fig6.add_subplot(111)

	allvals=[]			
	rectangle=[]

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
		plt.text(horizCount+half_inc,c+7,str(c),color='k',fontsize=fontsize,ha='center',va='center')
		plt.text(horizCount+half_inc,7,str(ConCount[i]),color='w',fontsize=fontsize,ha='center',va='center')
		#######


		plt.text(horizCount,bottomV,str(edges[i]),fontsize=15,rotation=90,ha='center',va='top')



		plt.plot([horizCount+.05,horizCount+.95],[bottomV,bottomV],linewidth=4,color='k')
		plt.plot([horizCount,horizCount],[0,bottomV],linewidth=3,color='k')


		horizCount+=30
		rowB=[]

	for r in rectangle:
		ax6.add_patch(r)
	plt.axis('equal')
	plt.axis('off')
	ax6.set_aspect('auto')
	plt.show()

def PlotSwarmSize():

	data["DeletionType"]=data["DeletionType"].str.replace('Capped','other')	
	##Capped is in place for very odd cases of deletions. Leave as others
	fig=plt.figure()
	ax=fig.add_subplot(111)
	sns.set(style="whitegrid",palette="muted")
	ax.set(yscale="log")
	ax=sns.swarmplot(x="production_centre",y="DeletionSize",hue="DeletionType",data=data,alpha=.6)
	axes=ax.axes
	axes.set_ylim(0,5100)	##cutting off some outliers
	plt.show()



if __name__ == "__main__":
	data=pd.read_pickle("crispr.clean.verified.aggregated.DF.pkl")
	FounderRate=data.g0_with_deletion_mutation/data.go_screened
	data["FounderRate"]=FounderRate

	#data.to_csv("./crispr.modified.csv",sep=";")##We already have comma-separated columns


	Sbins=[0,50,100,150,200,250,500,1000,2000,3000,5000]
	Fbins=list(np.arange(0,.9,0.1))


	bin_edges=map(lambda i: float(i),stats.mstats.mquantiles(data.DeletionSize,[np.arange(0,1,0.1)])[0])
	bin_edges[0]=0

	f_edges=map(lambda i: float(i),stats.mstats.mquantiles(data.FounderRate,[np.arange(0,1,0.1)])[0])
	f_edges[0]=0

	sb,Db=BinItems(data,Sbins,Fbins)
	Qb=BinFounders(data,Sbins)

	counts=[]
	edges=[]
	FounderR=[]
	ConCount=[]
	for k,j in sb.iteritems():
		FounderR.append(np.mean(j))
		counts.append(len(j))
		ConCount.append(len(Qb[k]))
		edges.append(Sbins[k])

	###########################################
	ScoresR=[]
	for k,j in Db.iteritems():
		if not len(j):	break
		ScoresR.append(np.mean(j))

	##########################################
	######## PLOTTING
	########
	PlotSizeEfficiencies()
	PlotSwarmSize()
