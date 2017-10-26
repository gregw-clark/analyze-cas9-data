#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn import svm
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import RFE
#from sklearn.feature_selection import SelectKBest
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


from sklearn.feature_extraction.text import CountVectorizer
from sklearn.svm import LinearSVC




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
	prod=1e9	##we have a good gaussian, but this scale gives is better 'looking' numbers. 
			##e.g. from 1.248e-9 to 1.24. All numbers end up as e-9, rescaled to >0 and < 10.
	for n in range(2,len(string)+1):
		if n >=18: argmin=4
		else:	argmin=min(n,3)
			
		cd=set(map(lambda q: "".join(q),list(window(string,n))))
		gg=list(map(lambda q: "".join(q),list(window(string,argmin))))
		
		prod*=len(cd)/float(len(gg))
#	if prod < 2:
#		print prod,string
#	elif prod > 7:
#		print "\t\t\t",prod,string
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


#df=pd.read_csv("cripsr_data_04_04_17.csv")

data=pd.read_pickle("Pared.DataFrame.pkl")

FounderRate=data.g0_with_deletion_mutation/data.go_screened
GLTRate=data.genotype_confirmed_f1s/data.g0_bred

data["GLTRate"]=GLTRate
data["FounderRate"]=FounderRate
olo=[]
for j in data.FounderRate:
	if j > 0.8:
		olo.append("Good")
	elif j <= 0.80 and  j > 0.25:
		olo.append("Fair")
	else:
		olo.append("Weak")
data["SuccessCode"]=olo
data=data[np.isfinite(data['FounderRate'])]


impBin=Imputer(strategy="most_frequent")
impFloat=Imputer(strategy="mean")

#binaryFeatures=["mrna_protein","mrna_nuclease_conc.values","secondary_allele","delivery_method","production_centre","grnas_per_gene"]
#binaryFeatures=['production_centre','primary_allele', 'secondary_allele', 'delivery_method', 'cas9_d10a', 'mrna_protein', 'mrna_nuclease_conc', 'grna_conc','DeletionType','SuccessCode']
binaryFeatures=['production_centre','delivery_method', 'cas9_d10a', 'mrna_protein','DeletionType','SuccessCode']

floatFeatures=['e_injected', 'embryos_survived', 'e_transferred', 'go_screened', 'total_count_of_mutagenized_g0_detected_by_screen', 'g0_with_nhej_mutation', 'g0_with_deletion_mutation', 'g0_bred', 'genotype_confirmed_f1s', 'DeletionSize', 'GLTRate', 'FounderRate']

####ols_to_norm = ['Age','Height']
#survey_data[cols_to_norm] = survey_data[cols_to_norm].apply(lambda x: (x - x.mean()) / (x.max() - x.min()))


Ff=data.filter(items=floatFeatures)

#Ff.dropna(inplace=True)
#Ff[normCols]=Ff[normCols].apply(lambda x: (x - x.mean()) / (x.max() - x.min()))

binary=data.filter(items=binaryFeatures)
binaryX=pd.get_dummies(binary)

#scaledX=pd.DataFrame(normalize(Ff,norm='l2'),columns=floatFeatures)

GI=Ff.index.values
binaryX=binaryX.reindex_axis(GI,axis=0)
lp=pd.concat([binaryX,Ff],axis=1)#,join_axes=[binaryX.index])
###SuccessCode_Fair', 'SuccessCode_Good', 'SuccessCode_Weak'
y=np.array(lp.SuccessCode_Good)
lp=lp.drop('SuccessCode_Good',axis=1)
lp=lp.drop('SuccessCode_Fair',axis=1)
lp=lp.drop('SuccessCode_Weak',axis=1)
lp=lp.drop('FounderRate',axis=1)
lp=lp.drop('GLTRate',axis=1)
print y


#print lp
#print list(lp)
#sys.exit()

#lp.dropna(inplace=True)
#print lp.shape
#y=np.array(data.)
#del lp['SuccessValue']
nms=list(lp)
X=[]
for index,row in lp.iterrows():
	addedRow=np.nan_to_num(np.array(row))
#	print max(addedRow),min(addedRow)
	X.append(addedRow)
X=np.vstack(X)


newNames=nms
newX=X
#####################################


feature_names=list(nms)


clf=ExtraTreesClassifier(n_estimators=1500,random_state=0)
clf=clf.fit(newX,y)
importances=clf.feature_importances_
std = np.std([tree.feature_importances_ for tree in clf.estimators_],axis=0)

#print X.shape
#model=SelectFromModel(clf,prefit=True)
#print model
#X_new=model.transform(X)
#print X_new.shape

indices = np.argsort(importances)[::-1]
#print std
for f in range(newX.shape[1]):
	if importances[indices[f]] > 0.01:#125:
		print f+1,newNames[indices[f]],importances[indices[f]]


#print X
#cv = CountVectorizer()
#cv.fit(X)
#print len(cv.vocabulary_)
#print cv.get_feature_names()
#X_train = cv.transform(data)
model = svm.SVC(kernel='linear', C=.8,gamma=1) 
# there is various option associated with it, like changing kernel, gamma and C value. Will discuss more # about it in next section.Train the model using the training sets and check score
model.fit(X, y)
print model.score(X, y)
coef=model.coef_.ravel()
#print coef

print len(coef),len(feature_names)
def getVals(coef,feature_names,top_features=20):
	#coef=classifier.coef_.ravel()
	top_positive_coefficients = np.argsort(coef)[-top_features:]
	top_negative_coefficients = np.argsort(coef)[:top_features]
	top_coefficients = np.hstack([top_negative_coefficients, top_positive_coefficients])
	plt.figure(figsize=(15, 5))
	#colors = ['red' if c < 0 else 'blue' for c in coef[top_coefficients]]
	colors = ['red' if c < 0 else 'blue' for c in coef[top_coefficients]]
	plt.bar(np.arange(2 * top_features), coef[top_coefficients], color=colors)
	feature_names = np.array(feature_names)
	plt.xticks(np.arange(1, 1 + 2 * top_features), feature_names[top_coefficients], rotation=60, ha='right')
	plt.show()

getVals(coef,feature_names)
	#print top_coefficients
	#print svm
#plot_coefficients(svm, cv.get_feature_names())


#model=SelectFromModel(clf,prefit=True)
#X_new=model.transform(X)
#print X_new.shape
#print X_new


#print ranking
# Plot pixel ranking
#plt.matshow(ranking, cmap=plt.cm.Blues)
#plt.colorbar()
#plt.title("Ranking of pixels with RFE")
#plt.show()




# The "accuracy" scoring is proportional to the number of correct
# classifications
#rfecv = RFECV(estimator=svc, step=1, cv=StratifiedKFold(2),
#              scoring='accuracy')
#rfecv.fit(X, y)

#print("Optimal number of features : %d" % rfecv.n_features_)

# Plot number of features VS. cross-validation scores
#plt.figure()
#plt.xlabel("Number of features selected")
#plt.ylabel("Cross validation score (nb of correct classifications)")
#plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
#plt.show()
