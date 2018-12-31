from glob import glob
import pandas as pd

def intersection(files, Chr):
	bedFiles = files
	if len(bedFiles) == 0:
		return []
	else:
		data = {"col1":[], "col2":[]}
		for file in bedFiles:
			f1 = open(file,"r")
			for line in f1:
				if Chr in line:
					aux = line[:-1].split("\t")
					data["col1"].append(int(aux[1]))
					data["col2"].append(int(aux[2]))
		df = pd.DataFrame(data=data)
		df = df.sort_values(by=["col1","col2"])
		finalData = {"pos1":[], "pos2":[], "rep":[]}
		previous = None	
		for index, row in df.iterrows():
			if previous != None:
				if previous[0] == row["col1"] and previous[1] == row["col2"]:
					previous[2] += 1
				else:
					finalData["pos1"].append(previous[0])
					finalData["pos2"].append(previous[1])
					finalData["rep"].append((previous[2]/len(bedFiles))*100)
					previous = [row["col1"],row["col2"],1]
			else:
				previous = [row["col1"],row["col2"],1]
	
		finalData["pos1"].append(previous[0])
		finalData["pos2"].append(previous[1])
		finalData["rep"].append((previous[2]/len(bedFiles))*100)
		dfFinal = pd.DataFrame(data = finalData)

		return dfFinal
		
def intersect(path_with_features, Chr):

	bedFiles = glob(path_with_features+"/*.bed")
	dfFinal = intersect(bedFiles, Chr)
	return dfFinal
