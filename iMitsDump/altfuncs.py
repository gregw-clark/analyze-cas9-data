#!/usr/bin.env python


import pandas as pd

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
def BinFounders(df,_S):
	SizeBins={0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}
	ScoreBins={0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}
	df=df[(df.status == "Founders Obtained") | (df.status == "Genotype confirmed")]
	for j,k in zip(df.DeletionSize,df.FounderRate):
		Sbin=whichBin(j,_S,"Score")
		SizeBins[Sbin].append(k)
	return SizeBins
