#!/usr/bin/env python
import vcf
import subprocess
# subprocess.Popen("vcftools --gzvcf /Users/beclose12/preQualDesktop/MSG_imp/SAU/chr11.dose.vcf.gz --snps HapListTest.txt --recode --recode-INFO-all --out sa_hap")
# subprocess.Popen("bgzip -c sa_hap.recode.vcf >sa_hap.vcf.gz")
# subprocess.Popen("tabix -p vcf sa_hap.vcf.gz")
## open file for writing Haps Classes
target = open('haploClasses_sa_pSnps.txt', 'w')
### read vcf.gz
vcf_reader = vcf.Reader(open('/Users/beclose12/preQualDesktop/HapClass/sa_pSnps.vcf.gz', 'r'))
## for the sake of testing use this
### SAU: rs3834466(5291563) = TT, rs7482144(5276169)=AA; rs28440105 (5269799)=CC; rs10128556 (5263683)= TT; rs968857 (5260458T) = TT  
 ## SAU: rs3834466 = GT (Hex), rs7482144=AA (allele2); rs28440105=CC (allele2); rs10128556= TT (allele1); rs968857 = TT  
 ## SEN: rs3834466 = G (Fam), rs7482144=AA (allele2); rs28440105=CC (allele2); rs10128556= TT (allele1); rs968857 = TT  
 ## BEN: rs3834466 = G (Fam), rs7482144=GG (allele1); rs28440105=CC (allele2); rs10128556= CC (allele1); rs968857 = TT  
 ## CAR: rs3834466 = G (Fam), rs7482144=GG (allele1); rs28440105=CC (allele2); rs10128556= CC (allele2); rs968857 = CC  
 ## CAM:   rs3834466 = G (Fam), rs7482144=GG (allele1); rs28440105=AA (allele1); rs10128556= CC (allele1); rs968857 = TT  
SAU_HET = {'rs3834466' : 'GT','rs28440105' : 'C', 'rs10128556' : 'T' , 'rs968857' : 'T'} 
SEN_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'T' , 'rs968857' : 'T'} 
BEN_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'C' , 'rs968857' : 'T'}
CAR_HET = {'rs3834466' : 'G','rs28440105' : 'C', 'rs10128556' : 'C' , 'rs968857' : 'C'}
CAM_HET =  {'rs3834466' : 'G','rs28440105' : 'A', 'rs10128556' : 'C' , 'rs968857' : 'T'} 
## get all samples
samples = vcf_reader.samples
## get all snp records
records = [record for record in vcf_reader]
## form hap classes dictionaries from the data
## for the sake of testing
for rec in records:
	print rec.ID
	print rec
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
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/SEN/SEN\/BEN/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/SAU/SAU\/BEN/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/CAR/CAR\/BEN/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/BEN\/CAM/CAM\/BEN/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/SAU\/SEN/SEN\/SAU/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/SEN\/CAR/CAR\/SEN/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/SEN\/CAM/CAM\/SEN/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/CAR\/SAU/SAU\/CAR/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/CAM\/SAU/SAU\/CAM/g", 'haploClasses_sa_pSnps.txt'])
sub = subprocess.call(['sed', '-i', '-e', r"s/CAM\/CAR/CAR\/CAM/g", 'haploClasses_sa_pSnps.txt'])