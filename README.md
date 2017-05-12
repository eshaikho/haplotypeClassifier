# haplotypeClassifier
Phased SNP-based Classification of Sickle Cell Anemia HBB Haplotypes.

This is a SNP-based Method to classify sickle cell haplotypes based on 4 SNPs.

## Getting Started

This scritp requires an input as vcf file of phased-imputed GWAS data.The phased GWAS data allow assigning of SNPs to paternal and maternal chromosomes which facilitates the classification procedure. Imputation and data preperation can be done using different tools and resources. We used Michigan Imputation server to impute data used in this project. See the link for more information.
([Michigan Imputation Server](https://imputationserver.sph.umich.edu/start.html#!pages/help)).
### Prerequisite

1. PyVCF ([PyVCF](https://github.com/jamescasbon/PyVCF)).

```
pip install pyvcf
```
          
2. Pysam ([Pysam](https://github.com/pysam-developers/pysam)).


```
pip install pysam
```
          
### Installing

Download the hapClassifier script and save it in the directory that contains the file you need to classify and run the following example commands:


**Run the classifier**
```
python hapClassifier.py input_file(bgzipped and tabix indexed file) output_file(the name of the file) 
```

## Author

* **Elmutaz Shaikho** - *Phased SNP-based Classification of Sickle Cell Anemia HBB Haplotypes* - [haplotypeClassifier](https://github.com/eshaikho/haplotypeClassifier)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* PyVCF developer 
* pysam developer
