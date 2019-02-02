from Bio import SeqIO
from numpy import zeros

def  generateDNAFragments(genome, fragmentSize, currentCHR):
	for record in SeqIO.parse(genome, "fasta"):
		if record.id == currentCHR:
			fragmentsChr=len(record.seq)/fragmentSize
			if int(fragmentsChr)<fragmentsChr:
				#for the current chr, if division is a float, the number of fragments will be the integer part + 1
				fragmentsChr = int(fragmentsChr)+1
			return zeros(int(fragmentsChr))

def fullVector(genome , currentCHR):
	for record in SeqIO.parse(genome, "fasta"):
		if record.id == currentCHR:
			return zeros(len(record.seq))
if __name__ == "__main__":
	print("Only available as library.\n\nexiting")
