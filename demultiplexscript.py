#!/usr/bin/env python3

from sys import version
import gzip
#print(version)
import argparse
import matplotlib

def get_args():
	parser = argparse.ArgumentParser(description="organize table of renamed file and if it's index or biological file")
	parser.add_argument("-f1", "--file1", help="biological forward reads", required=True)
	parser.add_argument("-f2", "--file2", help="index forward reads", required=True)
	parser.add_argument("-f3", "--file3", help="index reverse reads", required=True)
	parser.add_argument("-f4", "--file4", help="biological reverse reads", required=True)
	#parser.add_argument("-l", "--length", type=int, help="length of read; the index or biological read length", required=True)
	return parser.parse_args()
	
args = get_args()
f1 = args.file1
f2 = args.file2
f3 = args.file3
f4 = args.file4
#print("Using files:", f1, f2, f3, f4)

#l = args.length
file1 = [0, 0, 0, 0]
file2 = [0, 0, 0, 0]
file3 = [0, 0, 0, 0]
file4 = [0, 0, 0, 0]
#initializes four arrays for each of the four files with "0" placeholders; in each of the files the respective header, seqline, plus, and qscore line will be stored

barcodes = []

count = 0
undeter_ct = 0
hop_ct = 0
unhop_ct = 0


unhop_dict = {}

complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def rev_comp(DNAseq):
	'''Takes DNA molecule and returns the antisense complement'''
	complement = []
	for i in DNAseq:
		if i == "T":
			complement.append("A")
		if i == "A":
			complement.append("T")
		if i == "C":
			complement.append("G")
		if i == "G":
			complement.append("C")
		my_str = ''.join(complement)
	return my_str[::-1]

def convert_phred(letter):
	"""Converts a single character into a phred score"""
	return ord(letter)-33

def check_qual(seq):
	sum_qual = 0
	for char in seq:
		current_score = convert_phred(char) 
		sum_qual += current_score
	return sum_qual / len(seq)

fh = "/home/apulvino/bgmp/Bi622/DemultiplexProj/index_list.txt"

with open("index_list.txt", "r") as fh:
    for line in fh:
        line = line.strip()
        barcodes.append(line)
        barcodes.append(rev_comp(line))   
    #stores all of the known barcodes and their reverse complements in the list "barcodes"        
    #print(barcodes)  
   
    for barcode in barcodes:
        unhop_dict[barcode] = 0
	
	#stores all of the known barcodes and their reverse complements in the list "barcodes"		
	#for the counter in my dictionary (which is the value), let the dictionary key be the barcodes from the barcodes list to be stored in 
	#the variable unhop_ct (...the value); if the barcodes in the unhopped dictionary add one to the counter if not, append it to the dictionary

with gzip.open(f1, "r") as fh1, gzip.open(f2, "r") as fh2, gzip.open(f3, "r") as fh3, gzip.open(f4, "r") as fh4:
	
	while True:
		
		file1[0] = fh1.readline().decode('UTF-8').strip()
		file1[1] = fh1.readline().decode('UTF-8').strip()
		file1[2] = fh1.readline().decode('UTF-8').strip()
		file1[3] = fh1.readline().decode('UTF-8').strip()
		if file1[0] == "":
			break
		file2[0] = fh2.readline().decode('UTF-8').strip()
		file2[1] = fh2.readline().decode('UTF-8').strip()
		file2[2] = fh2.readline().decode('UTF-8').strip()
		file2[3] = fh2.readline().decode('UTF-8').strip()
		
		file3[0] = fh3.readline().decode('UTF-8').strip()
		file3[1] = fh3.readline().decode('UTF-8').strip()
		file3[2] = fh3.readline().decode('UTF-8').strip()
		file3[3] = fh3.readline().decode('UTF-8').strip()
		
		file4[0] = fh4.readline().decode('UTF-8').strip()
		file4[1] = fh4.readline().decode('UTF-8').strip()
		file4[2] = fh4.readline().decode('UTF-8').strip()
		file4[3] = fh4.readline().decode('UTF-8').strip()
	
		undetermined_fw = open("undetermined_fw.fq", "a")
		undetermined_rv = open("undetermined_rv.fq", "a")
		hopped_fw = open("hopped_fw.fq", "a")
		hopped_rv = open("hopped_rv.fq", "a")
#reads in all the gzipped files, strips the newline characters and then open each of the hopped and undetermined files "appendably" ("a") so you can add stuff
		
		#print(file2[1])
		#print(file3[1])
		#print(check_qual(file1[1]))
		#print(check_qual(file2[1]))
		#print(check_qual(file3[1]))
		#print(check_qual(file4[1]))
		#print(file3[1], file2[1], rev_comp(file3[1]), rev_comp(file3[1]) == file2[1])
		
		
		if 'N' in file2[1] or 'N' in file3[1] or file2[1] not in barcodes or file3[1] not in barcodes or check_qual(file2[1])<26 or check_qual(file3[1])<26:
			undetermined_fw.write(f"{file1[0]} {file2[1]}-{file3[1]}")
			undetermined_fw.write('\n')
			undetermined_fw.write(file1[1])
			undetermined_fw.write('\n')
			undetermined_fw.write(file1[2])
			undetermined_fw.write('\n')
			undetermined_fw.write(file1[3])
			undetermined_fw.write('\n')
			undetermined_rv.write(f"{file4[0]} {file2[1]}-{file3[1]}")
			undetermined_rv.write('\n')
			undetermined_rv.write(file4[1])
			undetermined_rv.write('\n')
			undetermined_rv.write(file4[2])
			undetermined_rv.write('\n')
			undetermined_rv.write(file4[3])
			undetermined_rv.write('\n')
			undeter_ct += 1
			count += 1
			#if there are Ns in the barcodes files, or if the barcodes don't match the ones we know were used in the prep, or if they 
			# have dumpy quality scores (less than 26, around 90% good) then we shoved the records in undetermined files with the barcodes 
			# that are there appended to the header line 
			# also increments for the undetermined reads counter
			# increments for the overall counter
			
		elif rev_comp(file3[1]) == file2[1]:
			#print(unhop_ct, "how")
			unhopped_fw = open(f"{file2[1]}-{file3[1]}_unhopped_fw.fq", "a")
			unhopped_rv = open(f"{file2[1]}-{file3[1]}_unhopped_rv.fq", "a")
			unhopped_fw.write(f"{file1[0]} {file2[1]}-{file3[1]}")
			unhopped_fw.write('\n' + file1[1] + '\n' + file1[2] + '\n' + file1[3] + '\n')
			unhopped_rv.write(f"{file1[0]} {file2[1]}-{file3[1]}")
			unhopped_rv.write('\n' + file4[1] + '\n' + file4[2] + '\n' + file4[3] + '\n')
			## if the reverse complement of the reverse complement barcode is the same as the sense barcode then open the
			## unhoppped/matched files and write the records out to them
			unhop_dict[file2[1]] += 1
			unhop_dict[file3[1]] += 1
			unhop_ct += 1
			count += 1
			#adds 1 when barcodes that match are iterated over in the dictionary for sense and antisense indexes
			# increments the counter for both unhopped and an overall count (for all reads in all files)
			
			unhopped_fw.close()
			unhopped_rv.close()	
			##then close those files
			
		else:
			#print(hop_ct, "when")
			hopped_fw.write(f"{file1[0]} {file2[1]}-{file3[1]}")
			hopped_fw.write('\n' + file1[1] + '\n' + file1[2] + '\n' + file1[3] + '\n')
			hopped_rv.write(f"{file1[0]} {file2[1]}-{file3[1]}")
			hopped_rv.write('\n' + file4[1] + '\n' + file4[2] + '\n' + file4[3] + '\n')
			hop_ct += 1
			count += 1
			## writes out all the hopped records with their hopped indexes stored in the header line
			
			
			## if the above conditionals are not met then we call these records index hopped and we write them out to the hopped file and increment our
			## counter (initialized at the very top) by one
			hopped_fw.close()
			hopped_rv.close()
			undetermined_fw.close()
			undetermined_rv.close()
			## closes the undetermined and hopped files
		
		#if count == 100000:
		#	break
		##checks the first 10,0000 lines to make sure everything is counting smoothly
		
for key, value in unhop_dict.items():
	print(key, value, count, (value/count)*100)
print("Percentage of read pairs with properly matched indexes:", unhop_ct/count*100)
print("Percentage of read pairs with improperly matched indexes:", hop_ct/count*100)
print("Percentage of read pairs with unknown indexes:", undeter_ct/count*100)
## prints out the value of the counters divided by 100 to give us a picture of the fraction of reads with matched, improperly matched/hopped, or unknown indexes