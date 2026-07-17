# Test Hardy-Weinberg Equilibrium (HWE) on the X chromosome
The robust allele-based regression framework for testing Hardy-Weinberg equilibrium on X chromosome

## Usage

Download [RA_Xchr_HWE.R](Xchr_HWE/RA_Xchr_HWE.R) to a local folder.


```r
source('RA_Xchr_HWE.R')
```


## Non-pseudoautosomal Region (NPR)

In Xchr NPR, most females are diploid while most males are hemizygous. Suppose the variant of interest is bi-allelic and A is the reference allele. Suppose we have a homogeneous population, female genotype counts are denoted as $f_{AA}$, $f_{AB}$ and $f_{BB}$; while male genotype counts are denoted as $m_A$ and $m_B$. 

### Input to the 'RA_Xchr' function. 

### Joint test of Hardy-Weinberg disequilibrium and sex-difference in minor allele frequency (sdMAF)

```r
RA_Xchr(gF, gM, snp_type='NPR', sdMAF, use_male)
```


```r
RA_Xchr(gF, gM, snp_type='NPR', sdMAF, use_male)
```

### X chromosome PAR SNPs

```r
RA_Xchr(gF, gM, snp_type='PAR', sdMAF, use_male)

```
