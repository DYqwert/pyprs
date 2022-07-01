# pyprs
Given the weights of different loci, the polygenic risk score of existing yeast was calculated
Manual
========

PRSCalculater_hybrid_SNP.py

A parser for calculating the individual polygenic risk score of F1 hybrids that accepts
input from the user

optional arguments:
  -h, --help           show this help message and exit
  --anno ANNO          GCF file for finding gene name
  --riskloci RISKLOCI  GWAS following sig loci
  --snpdir SNPDIR      a directory of snpfiles
  --output OUTPUT      filename for output

PRSCalculater_SNP.py

A parser for calculating the individual polygenic risk score only based on SNPs of isolates  that accepts
input from the user

optional arguments:
  -h, --help           show this help message and exit
  --anno ANNO          GCF file for finding gene name
  --cnv CNV            CNV file of individuals
  --riskloci RISKLOCI  GWAS following sig loci
  --snpdir SNPDIR      a directory of snpfiles
  --output OUTPUT      filename for output

PRSCalculater_2.py

A parser for calculating the individual polygenic risk score only based on SNPs and CNVs of isolates that accepts
input from the user

optional arguments:
  -h, --help           show this help message and exit
  --anno ANNO          GCF file for finding gene name
  --cnv CNV            CNV file of individuals
  --riskloci RISKLOCI  GWAS following sig loci
  --snpdir SNPDIR      a directory of snpfiles
  --output OUTPUT      filename for output
