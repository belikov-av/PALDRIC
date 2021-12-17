import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from utils.utils import  sumVect, cancers, genders, rows

def step17(algo):
	print ("step 17 start")
	input_file = algo[:-4] + '/data/'+algo[:-4] +'_patients.tsv'
	algo_name = algo[:-4]

	from textwrap import wrap
	try:
		os.mkdir(algo_name)
	except:
		pass
	try:
		os.mkdir(algo_name+"/cumulative histograms")
	except:
		pass

	header = "Tumor stage"+"\t"+"\t".join(rows)+"\n"

	stages1 = ['I']#'I or II NOS'
	stages2 = ['II']
	stages3 = ['III']
	stages4 = ['IV']

	total_events = {}
	male_events ={}
	female_events = {}


	with open(input_file,'r') as inpt:
		inpt.readline().split("\t")
		
		for line in inpt:
			s = line.split("\t")
			if s[4] in stages1:
				key = '1'
			if s[4] in stages2:
				key = '2'
			if s[4] in stages3:
				 key = '3'
			if s[4] in stages4:
				key = '4'
			if key not in total_events.keys():
				total_events[key] = [0 for i in range(11)]
				total_events[key].append(set())
			total_events[key][-1].add(s[0]) #add patient
			for i in range(11):
				total_events[key][i] += int(s[5+i])
			if s[2] == "MALE":
				if key not in male_events.keys():
					male_events[key] = [0 for i in range(11)]
					male_events[key].append(set())
				male_events[key][-1].add(s[0]) #add patient
				for i in range(11):
					male_events[key][i] += int(s[5+i])
			###females
			if s[2] == "FEMALE":
				if key not in female_events.keys():
					female_events[key] = [0 for i in range(11)]
					female_events[key].append(set())
				female_events[key][-1].add(s[0]) #add patient
				for i in range(11):
					female_events[key][i] += int(s[5+i])

	keys = ['1','2','3','4']
	with open(algo_name+'/data/'+algo_name+"_distribution_stages.tsv",'w') as out:
		out.write(header)
		stages_plot = {}
		for c in keys:
			out.write(c+"\t")
			num = len(total_events[c][-1])
			a = total_events[c][:-1]
			a = list(map(lambda x:x/num,a))
			stages_plot[c] = a
			line = "\t".join(list(map(lambda x:str(x),a)))+"\n"
			out.write(line)

	x = keys
	data = stages_plot
	fig, ax = plt.subplots()
	title = 'Driver event distribution by cancer stage'
	for i in x:
		ax.bar(i, data[i][0],color = 'lightcoral')
		bott = data[i][0]
		ax.bar(i, data[i][1], bottom = bott,color = 'red')
		bott+= data[i][1]
		ax.bar(i, data[i][2], bottom = bott,color = 'maroon')
		bott+= data[i][2]
		ax.bar(i, data[i][3], bottom = bott,color = 'cyan')
		bott+= data[i][3]
		ax.bar(i, data[i][4], bottom = bott,color = 'c')
		bott+= data[i][4]
		ax.bar(i, data[i][5], bottom = bott,color = 'teal')
		bott+= data[i][5]
		ax.bar(i, data[i][6], bottom = bott,color = 'orange')
		bott+= data[i][6]
		ax.bar(i, data[i][7], bottom = bott,color = 'darkorange')
		bott+= data[i][7]
		ax.bar(i, data[i][8], bottom = bott,color = 'violet')
		bott+= data[i][8]
		ax.bar(i, data[i][9], bottom = bott,color = 'purple')
	ax.set_ylabel('Average number driver of events', fontdict = {'fontsize':30})
	ax.set_xlabel('Cancer stage', fontdict = {'fontsize':30})
	ax.set_title(title, pad = 18,fontdict = {'fontsize':30})
	plt.xticks(fontsize = 19)
	plt.yticks(fontsize = 19)
	fig.set_figwidth(10)
	fig.set_figheight(16)
	fig.set_facecolor('floralwhite')
	ax.set_facecolor('seashell')
	plt.savefig(algo_name+"/"+"/cumulative histograms/"+algo_name+"_"+"distribution_stages"+".pdf",dpi = 600)

	with open(algo_name+'/data/'+algo_name+"_distribution_stages_males.tsv",'w') as out:
		out.write(header)
		males_plot = {}
		for c in keys:
			out.write(c+"\t")
			num = len(male_events[c][-1])
			a = male_events[c][:-1]
			a = list(map(lambda x:x/num,a))
			males_plot[c] = a
			line = "\t".join(list(map(lambda x:str(x),a)))+"\n"
			out.write(line)
	x = keys
	data = males_plot
	fig, ax = plt.subplots()
	title = 'Driver event distribution by cancer stage in males'
	for i in x:
		ax.bar(i, data[i][0],color = 'lightcoral')
		bott = data[i][0]
		ax.bar(i, data[i][1], bottom = bott,color = 'red')
		bott+= data[i][1]
		ax.bar(i, data[i][2], bottom = bott,color = 'maroon')
		bott+= data[i][2]
		ax.bar(i, data[i][3], bottom = bott,color = 'cyan')
		bott+= data[i][3]
		ax.bar(i, data[i][4], bottom = bott,color = 'c')
		bott+= data[i][4]
		ax.bar(i, data[i][5], bottom = bott,color = 'teal')
		bott+= data[i][5]
		ax.bar(i, data[i][6], bottom = bott,color = 'orange')
		bott+= data[i][6]
		ax.bar(i, data[i][7], bottom = bott,color = 'darkorange')
		bott+= data[i][7]
		ax.bar(i, data[i][8], bottom = bott,color = 'violet')
		bott+= data[i][8]
		ax.bar(i, data[i][9], bottom = bott,color = 'purple')

	ax.set_ylabel('Average number of driver events', fontdict = {'fontsize':30})
	ax.set_xlabel('Cancer stage', fontdict = {'fontsize':30})
	ax.set_title("\n".join(wrap(title,35)), pad = 18,fontdict = {'fontsize':30})
	plt.xticks(fontsize = 19)
	plt.yticks(fontsize = 19)
	fig.set_figwidth(10)
	fig.set_figheight(16)
	fig.set_facecolor('floralwhite')
	ax.set_facecolor('seashell')
	plt.savefig(algo_name+"/"+"/cumulative histograms/"+algo_name+"_"+"distribution_stages_males"+".pdf",dpi = 600)

	with open(algo_name+'/data/'+algo_name+"_distribution_stages_females.tsv",'w') as out:
		females_plot = {}
		out.write(header)
		for c in keys:
			out.write(c+"\t")
			num = len(female_events[c][-1])
			a = female_events[c][:-1]
			a = list(map(lambda x:x/num,a))
			females_plot[c] = a
			line = "\t".join(list(map(lambda x:str(x),a)))+"\n"
			out.write(line)
	x = keys
	data = females_plot
	fig, ax = plt.subplots()
	title = 'Driver event distribution by cancer stage in females'
	for i in x:
		ax.bar(i, data[i][0],color = 'lightcoral')
		bott = data[i][0]
		ax.bar(i, data[i][1], bottom = bott,color = 'red')
		bott+= data[i][1]
		ax.bar(i, data[i][2], bottom = bott,color = 'maroon')
		bott+= data[i][2]
		ax.bar(i, data[i][3], bottom = bott,color = 'cyan')
		bott+= data[i][3]
		ax.bar(i, data[i][4], bottom = bott,color = 'c')
		bott+= data[i][4]
		ax.bar(i, data[i][5], bottom = bott,color = 'teal')
		bott+= data[i][5]
		ax.bar(i, data[i][6], bottom = bott,color = 'orange')
		bott+= data[i][6]
		ax.bar(i, data[i][7], bottom = bott,color = 'darkorange')
		bott+= data[i][7]
		ax.bar(i, data[i][8], bottom = bott,color = 'violet')
		bott+= data[i][8]
		ax.bar(i, data[i][9], bottom = bott,color = 'purple')
	ax.set_ylabel('Average number of driver events', fontdict = {'fontsize':30})
	ax.set_xlabel('Cancer stage', fontdict = {'fontsize':30})
	ax.set_title("\n".join(wrap(title,35)),pad = 18, fontdict = {'fontsize':30})
	plt.xticks(fontsize = 19)
	plt.yticks(fontsize = 19)
	fig.set_figwidth(10)   
	fig.set_figheight(16)  
	fig.set_facecolor('floralwhite')
	ax.set_facecolor('seashell')
	plt.savefig(algo_name+"/"+"/cumulative histograms/"+algo_name+"_"+"distribution_stages_females"+".pdf",dpi = 600)
	print ("step 17 start")