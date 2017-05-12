#!/usr/bin/env python
################################################################################
# hapClassifier.py is script to classify sickle cell haplotypes	
# Copyright (c) [2017] [Elmutaz Shaikho Elhaj Mohammed]
################################################################################
#!/usr/local/bin/
import vcf
import subprocess
import sys
import pysam
import os
## This script takes the input file as first argument and output file as second argument
## python hapClassifier_rs_pSnps_git.py input(bgzipped and tabix indexed file) output file(text file name)
## Check input and out put files
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
vcf_reader = vcf.Reader(filename=inputFile)	

## Get all samples
samples = vcf_reader.samples

## Define a list to fetch SNPs by position. coordinates are in the zero-based
## Coordinates Should be according to Human Genome Reference b37
snpList = [5291562, 5269798, 5263682, 5260457]

## get 4 SNPs records needed for classification
records=[record for i in snpList for record in vcf_reader.fetch('11', i,(i+1))]

## Check in records is empty
if not records:
	print "Error: Data may does not have the SNPs need for classification,\n \
	or SNPs locations and/or RSID do not match Humanan genome Reference b37 dbSNP build 141"
	exit(1)

## Check if the data is phased or not, if not EXIT and print error messege
call = record.genotype(samples[0])
if not call.phased:
	print "Error: Data is not phased"
	exit(1)
if "rs" in records[1].ID:
	AI_HET = {'rs3834466' : 'GT','rs28440105': 'C', 'rs10128556' : 'T' , 'rs968857' : 'T'} 
	SEN_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'T' , 'rs968857' : 'T'} 
	BEN_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'C' , 'rs968857' : 'T'}
	CAR_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'C' , 'rs968857' : 'C'}
	CAM_HET = {'rs3834466' : 'G','rs28440105' : 'A', 'rs10128556' : 'C' , 'rs968857' : 'T'}  
	print 'success'
else:
	AI_HET = {'11:5291563' : 'GT','11:5269799': 'C', '11:5263683' : 'T' , '11:5260458' : 'T'} 
	SEN_HET = {'11:5291563' : 'G','11:5269799' : 'C', '11:5263683' : 'T' , '11:5260458' : 'T'} 
	BEN_HET = {'11:5291563' : 'G','11:5269799' : 'C', '11:5263683' : 'C' , '11:5260458' : 'T'}
	CAR_HET = {'11:5291563' : 'G','11:5269799' : 'C', '11:5263683' : 'C' , '11:5260458' : 'C'}
	CAM_HET = {'11:5291563' : 'G','11:5269799' : 'A', '11:5263683' : 'C' , '11:5260458' : 'T'}

## Form hap classes dictionaries from the data		
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
		elif Hapclass_M == Hapclass_F == AI_HET:
			target.write( '\t'.join( [smpl, 'AI HOMO'] ) + '\n' )
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
		elif Hapclass_F == AI_HET:
			classF = 'AI'		
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
		elif Hapclass_M == AI_HET:
			classM = 'AI'		
		elif Hapclass_M == CAR_HET:
			classM = 'CAR'		
		elif Hapclass_M == CAM_HET:
			classM = 'CAM'		
		else:
			classM = 'UNK'
	if Hapclass_M != Hapclass_F:
		target.write( smpl +'\t' + classF + '/' + classM  + '\n' )					
target.close()

## Change hetro classification eg BEN/CAM and CAM/BEN to be only one of them BEN/CAM
sub = subprocess.call(['sed', '-i.bak', r"s/BEN\/SEN/SEN\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/BEN\/AI/AI\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/BEN\/CAR/CAR\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/BEN\/CAM/CAM\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/AI\/SEN/SEN\/AI/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/SEN\/CAR/CAR\/SEN/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/SEN\/CAM/CAM\/SEN/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/CAR\/AI/AI\/CAR/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/CAM\/AI/AI\/CAM/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/CAM\/CAR/CAR\/CAM/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/BEN\/UNK/UNK\/BEN/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/SEN\/UNK/UNK\/SEN/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/CAM\/UNK/UNK\/CAM/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/CAR\/UNK/UNK\/CAR/g", outputFile])
sub = subprocess.call(['sed', '-i.bak', r"s/AI\/UNK/UNK\/AI/g", outputFile])

## Remove backup file
os.remove(outputFile+r".bak")
