#importing libraries

import argparse   # arguments parser
from glob import glob
import sys
import os
from Bio import SeqIO
import tempfile
import random

# importing my own libraries
sys.path.append(str(os.path.realpath(__file__))[:-8]+"/libs")
from getFeatures import getFeatures
from getFeatures import getTFBSac

if __name__=="__main__":

	#input
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--genome", help="Path to genome file (in fasta format)", required = True)
	parser.add_argument("-f","--fragment_size", help="number of bases to fragment DNA sequence. Default: 50", type = float, default=50)
	parser.add_argument("-p","--percentage", help="For some experiments you may need to merge them, so to do it you need to incluye a percentage of files that contain the feature. Default: 100", type = float, default = 100)
	parser.add_argument("-o","--output", help="folder where results will be placed. Default: ./", default = "./")
	parser.add_argument("-t","--tmp", help="folder for temporal files. Default: "+tempfile.gettempdir(), default = tempfile.gettempdir())

	####################################
	##
	## DNase-seq arguments
	##
	####################################
	parser.add_argument("-ds","--DNase",help="Option to generate features for DNAse.", action="store_true")
	parser.add_argument("-dsf","--DNase_folder", help="Folder where dnaseq results for DHS (in bed file) are located.")
	parser.add_argument("-dsfm","--DNase_filling_mode", help="If you wish your data in  binary (is or not present) or in number (number of the characteristic). Options are \"binary\", \"normalize\"  or \"decimal\". Default: binary", default="binary")

	####################################
	##
	## DNA methylation arguments
	##
	####################################
	parser.add_argument("-dm","--dna_methylation",help="Option to generate features for dna methylation.", action="store_true")
	parser.add_argument("-dname","--methylation_Folder", help="Folder where dnaseq results for DHS (in bed file) are located.")
	parser.add_argument("-dmfm","--DNA_meth_filling_mode", help="If you wish your data in  binary (is or not present) or in number (number of the characteristic). Options are \"binary\", \"normalize\", \"percentage\"  or \"decimal\". Default: binary", default="binary")

	####################################
	##
	## Histone marks arguments
	##
	####################################
	parser.add_argument("-hm","--histonemarks",help="Option to generate features for histone marks.", action="store_true")
	parser.add_argument("-hmf","--histoneMarksFolder", help="Folder where histone marks ChIP-seq results (in bed file) are located.")
	parser.add_argument("-hmfm","--histoneMarks_filling_mode", help="If you wish your data in  binary (is or not present) or in number (number of the characteristic). Options are \"binary\", \"normalize\"  or \"decimal\". Default: binary", default="binary")
	#parser.add_argument("-hmbt","--histoneMarks_by_type", help="If you wish to get all histone marks in one vector or different vector of informations for histone mark. Options are together  or unique. Default: unique", default="unique")

	####################################
	##
	## TFBSac arguments
	##
	####################################
	parser.add_argument("-tf","--tf_folder", help="Folder where tf subfolders with ChIP-seq results (in bed file) are located.", required = True)
	parser.add_argument("-tft","--tf_type", help="Type of tf result. Values can be \"optimal\" or \"conservative\". Default: optimal", default="optimal")
	parser.add_argument("-ff", "--fimo_folder", help="Path to fimo folder. this folder have a subfolder for each chromosome, and in that subfolder is placed the fimo.txt result", required = True)
	parser.add_argument("-ffi","--fimo_filter", help="Number equal or higher of motifs in the database to consider it as a motif. Default: 20", default=20, type=int)

	#parsing arguments
	args = parser.parse_args()
	tempFilenames = [] # to save temporal data
	abc = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
	ID ="" # a random identifier to save data
	for i in range (10):
		ID += random.choice("abcde")


	if args.fragment_size <= 0:
		print ("The number of kb to parse dna (--DNA_split option) must be higher than 0. Exiting")
		exit()
	#for each chr
	for record in SeqIO.parse(args.genome, "fasta"):
		
		currentCHR = "chr"+record.id
		currentCHR = "chr20" #delete or comment this line  for full dataset creation toghether with the exit line at the end of the script
		############################################
		##
		## Getting TF list with experimental data
		## so in a future step I can look if exist
		## the TFBM to compare if there is an
		## TBBSac
		##
		############################################
		tfList = []
		tfCheckpoint = open("./dat/TFCheckpoint_download_180515.txt","r")
		tfDict = {} # to save paths for each tf
		for line in tfCheckpoint:
			if "gene_symbol" not in line:
				tfList.append(line.split("\t")[0])
		aux = glob(args.tf_folder+"/*")
		for tfDir in aux:
			flag = False
			for files in glob(tfDir+"/*.bed"):
				if args.tf_type in files:
					flag = True
			if flag == True:
				tfName = tfDir.replace("eGFP-","").split("/")[-1]
				if tfName in tfList:
					if tfName not in tfDict:
						tfDict[tfName] = [glob(tfDir+"/*"+args.tf_type+"*.bed")[0]]
					else:
						tfDict[tfName].append(glob(tfDir+"/*"+args.tf_type+"*.bed")[0])
			del (flag)
		del (aux)	
		del(tfCheckpoint)
		del(tfList)
		
		###################################
		##
		## getting vectors
		##
		###################################

#		#for dnase
#		DNase = None
#		if args.DNase:
#			DNase = getFeatures(args.DNase_folder, args.DNase_filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR)		
#			tempFilenames.append(args.tmp+"/DNase_"+ID)
#			f = open(args.tmp+"/DNase_"+ID,"w")
#			f.write("DNase\n")
#			for i in DNase:
#				f.write(str(i)+"\n")
#			f.close()
#			del(DNase)
#
#		methylation = None
#		#for dna methylation
#		if args.dna_methylation:
#			methylation = getFeatures(args.methylation_Folder, args.DNA_meth_filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR)
#			tempFilenames.append(args.tmp+"/methylation_"+ID)
#			f = open(args.tmp+"/methylation_"+ID,"w")
#			f.write("methylation\n")
#			for i in methylation:
#				f.write(str(i)+"\n")
#			f.close()
#			del(methylation)
#		
#		#for histone marks
#		if args.histonemarks:
#			for HMfolder in glob(args.histoneMarksFolder+"/*"):
#				histone_mark = getFeatures(HMfolder, args.histoneMarks_filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR)
#				tempFilenames.append(args.tmp+"/"+HMfolder.split("/")[-1]+"_"+ID)
#				f = open(args.tmp+"/"+HMfolder.split("/")[-1]+"_"+ID,"w")
#				f.write(HMfolder.split("/")[-1]+"\n")
#				for i in histone_mark:
#					f.write(str(i)+"\n")
#				f.close()
#				del(histone_mark)

		#for active TFBS
		if args.tf_folder and args.fimo_folder:
			TFBSac = getTFBSac(tfDict, args.fimo_folder, args.fimo_filter, args.genome, currentCHR, args.fragment_size, args.percentage)
			tempFilenames.append(args.tmp+"/TFBSac_"+ID)
			f = open(args.tmp+"/TFBSac"+ID,"w")
			f.write("TFBSac\n")
			for i in TFBSac:
				f.write(str(i)+"\n")
			f.close()
			del(TFBSac)		
		
		
		
		exit()

	#	dictDNase = []
	#	hmVector = []
	#
	#	#for DNase-seq
	#	if args.DNase:
	#		dictDNase = getFeatures(args.DNase_folder, args.DNase_filling_mode, args.fragment_size, args.occupancy, args.percentage, args.genome)

			
		#for histone marks
	#	if args.histonemarks:
	#		hmVector = histoneFeaturedVector(args.histoneMarksFolder, args.histoneMarks_filling_mode, args.fragment_size, args.occupancy, args.percentage, args.histoneMarks_by_type, args.genome)
	
	#	TFBSAc = getTFBSac(args.fimo_folder, args.tf_folder, args.fragment_size, args.occupancy, args.percentage, args.fimo_filter, dictMethylation, dictDNase, hmVector, args.genome)


#python3 main.py -g genome.fasta -ds -dsf ./encode_K562/DNase-seq/ -dsfm percentage -dm -dname encode_K562/WGBS/ -dmfm percentage -hm -hmf ./encode_K562/HM/ -hmfm percentage -tf ./encode_K562/ChIP-seq/ -ff ./encode_K562/motifs/
