#!/usr/bin/env python

import argparse

def get_args():
	parser = argparse.ArgumentParser(description="organize table of renamed file and if it's index or biological file")
	parser.add_argument("-f", "--file", help="the index or biological file", required=True)
	parser.add_argument("-l", "--length", type=int, help="length of read; the index or biological read length", required=True)
	return parser.parse_args()
	
args = get_args()
f = args.file
l = args.length

import gzip
#gzip functionality for python to read zipped files
import matplotlib
import numpy as np

file = f.strip('/projects/bgmp/shared/2017_sequencing/')

	

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter)-33
#phred score converter
	#for file in fh:
def populate_array(file, length):
    mean_scores = np.zeros(shape=(int(length)))
    with gzip.open(file, "rt") as fh:
#file rt = read text... this is how you run the python gzip module on your files
        LN = 0
        for line in fh:
            LN +=1
            if LN % 4 == 0:
                j = 0	
                for char in line.strip():
                	current_score = convert_phred(char)
                	mean_scores[j] += current_score
                	j +=1
    return(mean_scores, LN)
    

mean_scores, LN = populate_array(f, l)    
print(mean_scores) 
#prints array of mean scores 


meanscores_array = []
meanscores_array.append(mean_scores/(LN/4))

print(meanscores_array)

import matplotlib.pyplot as plt

l = int(l)

x = list(range(l))
y = mean_scores
plt.bar(x, y) 

plt.title("mean scores of base pairs")

plt.xlabel("base pairs")
plt.ylabel("mean scores")

plt.savefig(file+".png")
            

