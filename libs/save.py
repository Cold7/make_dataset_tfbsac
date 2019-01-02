import pandas as pd

def saveDataset(files, output):
	dict = {}
	for file in files:
		f = open(file,"r")
		count = 0
		currentName = ""
		for line in f:
			if count == 0:
				currentName = line[:-1]
				dict[currentName] = []
				count = 1
			else:
				dict[currentName].append(line[:-1])
	dataFrame = pd.DataFrame(dict)
	dataFrame.to_csv(path_or_buf=output+"/output.tsv", sep='\t')
