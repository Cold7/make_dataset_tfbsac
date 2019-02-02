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

def getGeneFeatures(gtfFile, rnaFolder, filling_mode, fragmentSize, percentage, genome, currentCHR):
	gtf = open(gtfFile, "r")
	genes = {}
	#creating dict of genes
	for line in gtf:
		aux = line.split("\t")
		if aux[2] == "gene" and aux[0] == currentCHR:
			geneID = aux[8].split("\"")[1]
			genes[geneID] = {"init" : int(aux[3])-1, "end" : int(aux[4])-1, "count": 0} #-1 due it does not start from 0
	#adding count to count key
	rnaFiles = glob(rnaFolder+"/*.tsv")
	for file in rnaFiles:
		f = open(file,"r")
		for line in f:
			if "gene_id" not in line:
				aux = line[:-1].split("\t")
				if aux[0] in genes:
					genes[aux[0]]["count"] += float(aux[4])
	#getting the count average and doing full list of genes and count
	geneFullVector = fullVector(genome, currentCHR)
	countFullVector = fullVector(genome, currentCHR)
	for item in genes.items():
		genes[item[0]]["count"] /= len(rnaFiles)
		init = genes[item[0]]["init"]
		end = genes[item[0]]["end"]
		
		while init <= end:
			geneFullVector[init] = 1
			countFullVector[init] = genes[item[0]]["count"] 
			init += 1
	#doing fragments of n bases
	geneDNAfragment = generateDNAFragments(genome, fragmentSize, currentCHR)
	countsDNAfragment = generateDNAFragments(genome, fragmentSize, currentCHR)
	
	#filling DNAfragment
	for i in range(len(geneFullVector)):
		pos = int(i/fragmentSize)
		geneDNAfragment[pos] += geneFullVector[i]
		if countFullVector[i] > countsDNAfragment[pos]:
			countsDNAfragment[pos] = countFullVector[i]

	if filling_mode != "binary":
		if filling_mode == "normalize":
			for i in range(len(geneDNAfragment)):
				geneDNAfragment[i] = 1 - (1/(1 + geneDNAfragment[i]))
		if filling_mode == "percentage":
			for i in range(len(geneDNAfragment)):
				geneDNAfragment[i] = geneDNAfragment[i]/fragmentSize
	else:#binary case
		for i in range(len(geneDNAfragment)):
			if geneDNAfragment[i] > 0:
				geneDNAfragment[i] = 1	

	return [geneDNAfragment,countsDNAfragment]

def getTFBSac(tfDict, fimo_folder, fimo_filter, genome, currentCHR, fragmentSize, percentage, filling_mode, ensembl_gff, positive_threshold, negative_threshold):
	toReturn = []
	fimoDNAfragment = generateDNAFragments(genome, fragmentSize, currentCHR)
	ensemblDNAfragment = generateDNAFragments(genome, fragmentSize, currentCHR)
	dataTFBSac = fullVector(genome, currentCHR)	
	
	#############################################
	##
	## using fimo data
	##
	#############################################
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
	
	for tf in tfDict.items():
		tfName = tf[0]
		print("working with ", tfName)
		paths = tf[1]
		if tfName in motifFound: #if the current tf have a motif
			dataTFBS = np.zeros(len(dataTFBSac))
			print("\t", tfName, "have a motif")
			#getting tfbs predicted by fimo
			file = open(fimo_folder+"/"+currentCHR+"/fimo.txt")
			for line in file:
				if line[0] != "#" or line != "\n" or line != "":
					dataLine = line.split("\t")
					if dataLine[0].split("_")[0] == tfName:
						init = int(dataLine[3]) - 1
						end = int(dataLine[4]) - 1
						while init <= end:
							dataTFBS[init] = -1
							init += 1		
			#now for chip data
			features = intersection(paths, currentCHR)
			#marking position where the characteristic is present
			for index, row in features.iterrows():
				if percentage >= row['rep']:
					init = int(row['pos1']) - 1 # -1 because the count start from 0 and not for 1
					end = int(row['pos2']) - 1
					while init < end:
						dataTFBS[init] *= -1 #so motif used will be 1, not used -1 and no tfbs 0
						init += 1
					
			#definition of tfbs active or inactive due the presence or ausence of chip data into the current motif			
			for i in range(len(dataTFBS)):
				if dataTFBS[i] == -1:
					if dataTFBSac[i] != 1:
						dataTFBSac[i] = -1
				if dataTFBS[i] == 1:
					dataTFBSac[i] = 1
	
	#filling DNAfragment
	for i in range(len(dataTFBSac)):
		pos = abs(int(i/fragmentSize))
		fimoDNAfragment[pos] += dataTFBSac[i]

	if filling_mode != "binary":
		if filling_mode == "normalize":
			for i in range(len(fimoDNAfragment)):
				fimoDNAfragment[i] = 1 - (1/(1 + fimoDNAfragment[i]))
		if filling_mode == "percentage":
			for i in range(len(fimoDNAfragment)):
				fimoDNAfragment[i] = fimoDNAfragment[i]/fragmentSize
	else:#binary case
		for i in range(len(fimoDNAfragment)):
			if fimoDNAfragment[i] >= positive_threshold:
				if positive_threshold != 0:
					fimoDNAfragment[i] = 1
				else:
					if fimoDNAfragment[i] > 0:
						fimoDNAfragment[i] = 1
			if fimoDNAfragment[i] <= negative_threshold:
				if negative_threshold != 0:
					fimoDNAfragment[i] = -1
				else:
					if fimoDNAfragment[i] < 0:
						fimoDNAfragment[i] = -1
			
	toReturn.append(fimoDNAfragment)

	###########################################
	##
	## using ensembl data
	##
	###########################################
	dataTFBSac = fullVector(genome, currentCHR)	
	tfEnsemblIndex = {}
	regulatory_gff = open(ensembl_gff,"r")
	regDict = {}

	#getting TF sites from ensembl regulatory gff
	cont = 0
	for line in regulatory_gff:
		aux = line.split("\t")
		if "chr"+aux[0] == currentCHR and "feature_type=TF binding site" in line:
			init = int(aux[3]) -1
			end = int(aux[4]) - 1
			regDict[str(cont)] =  {"init": init, "end": end}
			cont += 1
			while init <= end:
				dataTFBSac[init] = -1
				init += 1
			
	regulatory_gff.close()
	del(regulatory_gff)

	#looping over all tfs with chipseq
	for tf in tfDict.items():
		tfName = tf[0]
		paths = tf[1]
		features = None
		#getting chip-seq positions
		try:
			features = intersection(tf[1], currentCHR)
			#marking position where the characteristic is present
			for index, row in features.iterrows():
				if row["rep"] == 100:
					init = int(row['pos1']) - 1 # -1 because the count start from 0 and not for 1
					end = int(row['pos2']) - 1
					if percentage >= row['rep']:
						while init < end:
							dataTFBSac[init] *= -1 #so those tfbs inactives will be -1, with no tfbs sites will be 0 and active tfbs will be 1
							init += 1
		except:
			pass

	#transforming the whole tfbs from ensemble into active or inactive (due some parts are not binding tfs)
	for item in  regDict.items():
		init = item[1]["init"]
		end = item[1]["end"]
		i = init
		ones = False
		while i < end:
			if dataTFBSac[i] == 1:
				ones = True
			i += 1
		if ones == True:
			i = init
			while i < end:
				dataTFBSac[i] = 1
				i += 1
				
	#filling DNAfragment
	for i in range(len(dataTFBSac)):
		pos = abs(int(i/fragmentSize))
		ensemblDNAfragment[pos] += dataTFBSac[i]


	if filling_mode != "binary":
		if filling_mode == "normalize":
			for i in range(len(ensemblDNAfragment)):
				ensemblDNAfragment[i] = 1 - (1/(1 + ensemblDNAfragment[i]))
		if filling_mode == "percentage":
			for i in range(len(ensemblDNAfragment)):
				ensemblDNAfragment[i] = ensemblDNAfragment[i]/fragmentSize
	else:#binary case
		for i in range(len(ensemblDNAfragment)):
			if ensemblDNAfragment[i] >= positive_threshold:
				if positive_threshold != 0:
					ensemblDNAfragment[i] = 1
				else:
					if ensemblDNAfragment[i] > 0:
						ensemblDNAfragment[i] = 1
					
			if ensemblDNAfragment[i] <= negative_threshold:
				if negative_threshold != 0: 
					ensemblDNAfragment[i] = -1
				else:
					if ensemblDNAfragment[i] <= 0:
						ensemblDNAfragment[i] = -1
		
	toReturn.append(ensemblDNAfragment)
	
	return toReturn
