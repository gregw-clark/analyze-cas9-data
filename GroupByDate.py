#!/usr/bin/env python

import string,re, sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data=open("AllImits.nan.csv",'r').readlines()
print data[0]
header=data[0].strip().split(";")
data.pop(0)
colnames=dict(zip(header,range(len(header))))
from collections import defaultdict,Counter


import matplotlib.dates as mdates

years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')



dated=defaultdict(list)
MIA=defaultdict(list)
#['', 'genotype_confirmed_date', 'total_count_of_mutagenized_g0_detected_by_screen', 'go_screened', 'production_centre', 'cas9_d10a', 'g0_with_nhej_mutation', 'e_injected', 'marker_symbol', 'mrna_protein', 'g0_with_deletion_mutation', 'grna_cleave_sites', 'grna_coordinates', 'g0_obtained_date', 'g0_bred', 'e_transferred', 'genotype_confirmed_f1s', 'grnas_per_gene', 'mrna_nuclease_conc', 'grna_sequences', 'secondary_allele', 'delivery_method', 'embryos_survived', 'primary_allele', 'grna_conc']
Ustat='Genotype confirmed','Micro-injection aborted'

symbols=[]
for m in range(len(data)):
	kl=data[m].strip().split(";")
	g0_date,genotype_confirmed,mi_date,centre,status,marker=kl[colnames['g0_obtained_date']],kl[colnames['genotype_confirmed_f1s']],kl[colnames['mi_date']],kl[colnames['production_centre']],kl[colnames['status']],kl[colnames['marker_symbol']]
#	print g0_date,genotype_confirmed,centre,
#Genotype confirmed

	if status  == 'Genotype confirmed': 
		#print g0_date
		month,year=g0_date.split("/")[1:]
		monthly=year+"_"+month
		dated[monthly].append(centre)
		year,month=mi_date.split("-")[:2]
		monthly=year+"_"+month
		MIA[monthly].append(centre)
		if year == "2017":
		#if mi_date.startswith("2017"):
			print mi_date	
	#else == 'Micro-injection aborted':
	else:
		year,month=mi_date.split("-")[:2]
		monthly=year+"_"+month
		MIA[monthly].append(centre)
	symbols.append(marker)
print list(set(symbols))
print len(set(symbols))
sys.exit()		
from datetime import datetime
dts=MIA.keys()
dts.sort()
tG=0
tA=0
Mi_attempts=[]
G_confirmed=[]
RD=[]
DATES=[]
for d in dts:
#	print d,len(datdatetime_object = datetime.strptime('Jun 1 2005  1:33PM', '%b %d %Y %I:%M%p')ed[d]),Counter(dated[d]),
	datetime_object = datetime.strptime(d, '%Y_%m')
	print datetime_object
	
	if not d.startswith("2013"):
		DATES.append(datetime_object)
		try:	
			dated[d]
			tG+=len(dated[d])
		except KeyError:
			pass	
		#print d,tG,tA
		tA+=len(MIA[d])
		print d,tG,tA
		
		RD.append(d)
		Mi_attempts.append(tA)
		G_confirmed.append(tG)
sys.exit()
fig, ax = plt.subplots()

#datemin = datetime.date(DATES[0], 1, 1)
#datemax = datetime.date(DATES[-1] + 1, 1, 1)
ax.set_xlim(DATES[0], DATES[-1])

ax.plot(DATES,Mi_attempts,linewidth=4,label="MI attempts")
ax.plot(DATES,G_confirmed,linewidth=4,label="Confirmed Genotype")
ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
#ax.format_ydata = Mi_attempts
ax.set_ylabel("# Injections")
legend = ax.legend(loc='upper left', shadow=True)
ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()



#plt.plot_ts(DATES,Mi_attempts,color='g',linewidth=2)
#plt.plot_ts(DATES,color='r',linewidth=2)
plt.show()		
