#!/usr/bin/env python

import string,re,sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
from cPickle import load,dump
import numpy as np
from collections import defaultdict
import math
from scipy.stats import ttest_ind,wilcoxon

blues=['#dae8f5', '#bad6ea', '#88bedc', '#539dcc', '#2a7ab9', '#0b559f']	# generated from sns.color_palette("Blues")
reds=['#fddbcb', '#fcaf93', '#fb8161', '#f44e38', '#d52221', '#a91016']		# generated from sns.color_palette("Reds")
blues.reverse()


fig6 = plt.figure(figsize=(800,60))
ax6 = fig6.add_subplot(111)
n=open('2data.txt','r').readlines()
counts=map(lambda l: int(l),n[0].split())
edges=map(lambda l: int(l),n[1].split())


allvals=[]			
rectangle=[]


lclCol=blues[2]
verticalCount=1
keepDates=[]
horizCount=0
for i,c in enumerate(counts):
	#print i,c	
#	rowA=[]
	#rectangle.append(patches.Rectangle((.3*horizCount*spacer,0.3*verticalCount),0.3*spacer,0.3,facecolor=lclCol,linewidth=0))
	rectangle.append(patches.Rectangle((horizCount,0),1,c,facecolor=lclCol,ec='w',linewidth=0.5))
	print horizCount,horizCount+1
	#plt.plot([(horizCount+.5)*.3*spacer,(horizCount+.5)*.3*spacer],[Vspot*.3+0.25,keepDates.count(k)*0.3+Vspot*0.3-0.25],color='k',linewidth=.8)
	plt.text(horizCount-.025,0,str(edges[i]),fontsize=10,rotation=90,ha='left',va='top')
	fontsize=max(8,c/17)
	plt.text(horizCount+.5,c+10,str(c),fontsize=fontsize,ha='center',va='center')
	horizCount+=1
	rowB=[]

for r in rectangle:
        ax6.add_patch(r)

def addyaxis(yprops,Mousetext):
	for k in range(len(Mousetext)):
		plt.text(0,.3*(k+1)+.15,Mousetext[k],yprops)
def addxaxis(xprops,IMPCtext):
	for k in range(len(IMPCtext)):
		plt.text(.3*(k+1)+.15,0.2,IMPCtext[k],xprops,rotation=280)

plt.axis('equal')
plt.axis('off')
plt.axis('tight')
plt.show()
#fig6.show()
#fig6.savefig('rect6.png', dpi=90, bbox_inches='tight')
