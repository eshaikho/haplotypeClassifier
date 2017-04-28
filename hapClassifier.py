#!/usr/bin/env python
import vcf
import subprocess
import sys
## This script takes the input file as first argument and output file as second argument
## python hapClassifier_rs_pSnps_git.py input(bgzipped and tabix indexed file) output file(text file name)
## This script takes the input file as first argument and output file as second argument
## python hapClassifier_rs_pSnps_git.py input(bgzipped and tabix indexed file) output file(text file name)

if len(sys.argv) < 3:
    print "Error: one or more arguement is/are missing"
    print "python hapClassifier_rs_pSnps_git.py input(bgzipped and tabix indexed file) output file(text file name)"
    exit(1)
else:
## Define input and output files	
	inputFile = sys.argv[1] 
	outputFile = sys.argv[2] 
	print 'Input file:'+ sys.argv[1]
	print 'Output file:'+ sys.argv[2]
 
## open file for writing Haps Classes
target = open(outputFile, 'w')
### read vcf.gz
vcf_reader = vcf.Reader(open(inputFile, 'r'))

## get all samples
samples = vcf_reader.samples
## get all snp records
records = [record for record in vcf_reader]

## Display the SNP records and print the ID to see whether they are in RSID fashion or position fashion
for rec in records:
	print rec.ID
	print rec
if "rs" in records[1].ID:
	SAU_HET = {'rs3834466' : 'GT','rs28440105': 'C', 'rs10128556' : 'T' , 'rs968857' : 'T'} 
	SEN_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'T' , 'rs968857' : 'T'} 
	BEN_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'C' , 'rs968857' : 'T'}
	CAR_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'C' , 'rs968857' : 'C'}
	CAM_HET = {'rs3834466' : 'G','rs28440105' : 'A', 'rs10128556' : 'C' , 'rs968857' : 'T'}  
	print 'success'
else:
	SAU_HET = {'11:5291563' : 'GT','11:5269799': 'C', '11:5263683' : 'T' , '11:5260458' : 'T'} 
	SEN_HET = {'11:5291563' : 'G','11:5269799' : 'C', '11:5263683' : 'T' , '11:5260458' : 'T'} 
	BEN_HET = {'11:5291563' : 'G','11:5269799' : 'C', '11:5263683' : 'C' , '11:5260458' : 'T'}
	CAR_HET = {'11:5291563' : 'G','11:5269799' : 'C', '11:5263683' : 'C' , '11:5260458' : 'C'}
	CAM_HET = {'11:5291563' : 'G','11:5269799' : 'A', '11:5263683' : 'C' , '11:5260458' : 'T'}

## form hap classes dictionaries from the data		
for smpl in samples:
	keys = keys =[records[i].ID for i in range(len(records))]
	values_M = [(records[i].genotype(smpl).gt_bases).split('|')[0] for i in range(len(records))]
	values_F = [(records[i].genotype(smpl).gt_bases).split('|')[1] for i in range(len(records))]
	Hapclass_M = dict(zip(keys, values_M))	
	Hapclass_F = dict(zip(keys, values_F))	
	classF=""
	classM=""
	if Hapclass_M == Hapclass_F:
		if Hapclass_M == Hapclass_F == BEN_HET:
			target.write( '\t'.join( [smpl, 'BEN HOMO'] ) + '\n' )
		elif Hapclass_M == Hapclass_F == SEN_HET:
			target.write( '\t'.join( [smpl, 'SEN HOMO'] ) + '\n' )
		elif Hapclass_M == Hapclass_F == SAU_HET:
			target.write( '\t'.join( [smpl, 'SAU HOMO'] ) + '\n' )
		elif Hapclass_M == Hapclass_F == CAR_HET:
			target.write( '\t'.join( [smpl, 'CAR HOMO'] ) + '\n' )
		elif Hapclass_M == Hapclass_F == CAM_HET:
			target.write( '\t'.join( [smpl, 'CAM HOMO'] ) + '\n' )
		else: 
			target.write( '\t'.join( [smpl, 'UNK HOMO'] ) + '\n' )
	###
	if Hapclass_F  and (Hapclass_M != Hapclass_F) :
		if Hapclass_F == BEN_HET:
			classF = 'BEN'		
		elif Hapclass_F == SEN_HET:
			classF = 'SEN'		
		elif Hapclass_F == SAU_HET:
			classF = 'SAU'		
		elif Hapclass_F == CAR_HET:
			classF = 'CAR'		
		elif Hapclass_F == CAM_HET:
			classF = 'CAM'		
		else:
			classF = 'UNK'		
	###
	if Hapclass_M and (Hapclass_M != Hapclass_F) :	
		if Hapclass_M == BEN_HET:
			classM = 'BEN'		
		elif Hapclass_M == SEN_HET:
			classM = 'SEN'		
		elif Hapclass_M == SAU_HET:
			classM = 'SAU'		
		elif Hapclass_M == CAR_HET:
			classM = 'CAR'		
		elif Hapclass_M == CAM_HET:
			classM = 'CAM'		
		else:
			classM = 'UNK'
	if Hapclass_M != Hapclass_F:
		target.write( smpl +'\t' + classF + '/' + classM  + '\n' )					
target.close()
## change hetro classification eg BEN/CAM and CAM/BEN to be only one of them BEN/CAM
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/SEN/SEN\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/SAU/SAU\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/CAR/CAR\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/CAM/CAM\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/SAU\/SEN/SEN\/SAU/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/SEN\/CAR/CAR\/SEN/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/SEN\/CAM/CAM\/SEN/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/CAR\/SAU/SAU\/CAR/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/CAM\/SAU/SAU\/CAM/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/CAM\/CAR/CAR\/CAM/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/UNK/UNK\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/SEN\/UNK/UNK\/SEN/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/CAM\/UNK/UNK\/CAM/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/CAR\/UNK/UNK\/CAR/g", outputFile])
sub = subprocess.call(['sed', '-i', '-e', r"s/SAU\/UNK/UNK\/SAU/g", outputFile])
