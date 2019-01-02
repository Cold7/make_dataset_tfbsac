from glob import glob
from intersectBed import intersect
from intersectBed import intersection
from fragment import generateDNAFragments
from fragment import fullVector
import numpy as np

def getFeatures(folder, filling_mode, fragmentSize, percentage, genome, currentCHR):
	data = fullVector(genome, currentCHR)
	features = intersect(folder, currentCHR)
	#marking position where the characteristic is present
	for index, row in features.iterrows():
		init = int(row['pos1']) - 1 # -1 because the count start from 0 and not for 1
		end = int(row['pos2']) - 1
		if percentage >= row['rep']:
			while init < end:
				data[init] = 1
				init += 1

	DNAfragment = generateDNAFragments(genome, fragmentSize, currentCHR)
	
	#filling DNAfragment
	for i in range(len(data)):
		pos = int(i/fragmentSize)
		DNAfragment[pos] += data[i]


	if filling_mode != "binary":
		if filling_mode == "normalize":
			for i in range(len(DNAfragment)):
				DNAfragment[i] = 1 - (1/(1 + DNAfragment[i]))
		if filling_mode == "percentage":
			for i in range(len(DNAfragment)):
				DNAfragment[i] = DNAfragment[i]/fragmentSize
	else:#binary case
		for i in range(len(DNAfragment)):
			if DNAfragment[i] > 0:
				DNAfragment[i] = 1
	return DNAfragment


def getTFBSac(tfDict, fimo_folder, fimo_filter, genome, currentCHR, fragmentSize, percentage, filling_mode):
	dataTFBSac = fullVector(genome, currentCHR)
	print("we need to change this in function of the TFBM DB used")
	motifFound = []
	fimoScanned = open(fimo_folder+"/"+currentCHR+"/fimo.txt","r")
	#we will look wich tf have tfbm
	for line in fimoScanned:
		if line[0] != "#":
			tfbm = line.split("\t")[0].split("_")[0]
			if tfbm not in motifFound:
				motifFound.append(tfbm)
	fimoScanned.close()
	del(fimoScanned)

	#looping over all tfs with chipseq
	for tf in tfDict.items():
		tfName = tf[0]
		paths = tf[1]
		if tfName in motifFound: #if the current tf have a motif
			#getting tfbs predicted by fimo
			dataTFBM = np.zeros(len(dataTFBSac))
			file = open(fimo_folder+"/"+currentCHR+"/fimo.txt")
			for line in file:
				if line[0] != "#" or line != "\n" or line != "":
					dataLine = line.split("\t")
					if dataLine[0].split("_")[0] == tfName:
						init = int(dataLine[3]) - 1
						end = int(dataLine[4]) - 1
						while init <= end:
							dataTFBM[init] = 1
							init += 1
			
			file.close()
			del(file)
			
			#getting chip-seq positions
			dataChIP = np.zeros(len(dataTFBSac))
			features = intersection(tf[1], currentCHR)
			#marking position where the characteristic is present
			for index, row in features.iterrows():
				if row["rep"] == 100:
					init = int(row['pos1']) - 1 # -1 because the count start from 0 and not for 1
					end = int(row['pos2']) - 1
					if percentage >= row['rep']:
						while init < end:
							dataChIP[init] = 1
							init += 1
				
			#definition of tfbs active or inactive due the presence or ausence of chip data into the current motif			
			for i in range(len(dataTFBM)):
				#if in dataTFBM and in dataChIP in the same position there exist a mark, there will we an TFBSactive
				if dataTFBM[i] == 1 and dataChIP[i] == 1:
					dataTFBSac[i] = 1
				if dataTFBM[i] == 1 and dataChIP[i] == 0:
					if dataTFBSac[i] != 1:
						dataTFBSac[i] == -1
			
	DNAfragment = generateDNAFragments(genome, fragmentSize, currentCHR)
	
	#filling DNAfragment
	for i in range(len(dataTFBSac)):
		pos = abs(int(i/fragmentSize))
		DNAfragment[pos] += dataTFBSac[i]


	if filling_mode != "binary":
		if filling_mode == "normalize":
			for i in range(len(DNAfragment)):
				DNAfragment[i] = 1 - (1/(1 + DNAfragment[i]))
		if filling_mode == "percentage":
			for i in range(len(DNAfragment)):
				DNAfragment[i] = DNAfragment[i]/fragmentSize
	else:#binary case
		for i in range(len(DNAfragment)):
			if DNAfragment[i] > 0:
				DNAfragment[i] = 1
			if DNAfragment[i] < 0:
				DNAfragment[i] = -1
	return DNAfragment
