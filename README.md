# haplotypeClassifier
Phased SNP-based Classification of Sickle Cell Anemia HBB Haplotypes.

This is a SNP-based Method to classify sickle cell haplotypes based on 4 SNPs.

## Getting Started

This scritp requires an input as vcf file of phased-imputed GWAS data.The phased GWAS data allow assigning of SNPs to paternal and maternal chromosomes which facilitates the classification procedure. Imputation and data preperation can be done using different tools and resources. We used Michigan Imputation server to impute data used in this project.See the link below for more information.
([Michigan Imputation Server](https://imputationserver.sph.umich.edu/start.html#!pages/help))
### Prerequisite

1. HTSlib that contains bgzip and tabix to zip and index the vcf files, respectively ([HTSlib](https://github.com/samtools/htslib)).
          
2. VCFtools to subset the SNPs from the imputed file ([VCFtools](https://vcftools.github.io/examples.html)).
          
3. PyVCF ([PyVCF](https://github.com/jamescasbon/PyVCF)).

### Installing

Download the hapClassifier script and save in the directory contains the file you need to classify and run the following example commands:

Prepare the input file
```
vcftools --gzvcf chr11.dose.vcf.gz --snps HapListTest.txt --recode --recode-INFO-all --out sa_hap
bgzip -c sa_hap.recode.vcf >sa_hap.vcf.gz
tabix -p vcf sa_hap.vcf.gz
```
Run the classifier
```
python hapClassifier.py input(bgzipped and tabix indexed file) output file(text file name)
```

## Authors

* **Elmutaz Shaikho** - *Initial work* - [hapClassifier](https://github.com/eshaikho/haplotypeClassifier)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* PyVCF developer 
* HTSliB developers
* VCFtools developers
