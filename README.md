# Test Hardy-Weinberg Equilibrium (HWE) on the X chromosome
The robust allele-based regression framework for testing Hardy-Weinberg equilibrium on X chromosome

## Usage

Download [RA_Xchr_HWE.R](Xchr_HWE/RA_Xchr_HWE.R) to a local folder.


```r
source('RA_Xchr_HWE.R')
```


## Non-pseudoautosomal Region (NPR)

In Xchr NPR, most females are diploid while most males are hemizygous. Suppose the variant of interest is bi-allelic and A is the reference allele. Suppose we have a homogeneous population, female genotype counts are denoted as $f_{AA}$, $f_{AB}$ and $f_{BB}$; while male genotype counts are denoted as $m_A$ and $m_B$. 

Use the `RA_Xchr` function to perform the test. Input arguments include
- `gF`: female genotype counts. A vector of length 3. `gF = c(fAA, fAB, fBB)`.
- `gM`: male genotype counts. A vector of length 2. `gM = c(mA, mB)`.
- `snp_type`: type of SNP. `snp_type = 'NPR' ` for Xchr NPR SNPs and `snp_type = 'PAR' ` for Xchr PAR SNPs.
- `sdMAF`: sdMAF status. Logical. If `snp_type = TRUE`, the reference allele frequency is estimated separately in females and males; if `snp_type = FALSE`, the reference allele frequency is estimated using a pooled sample of females and males. If the objective is to perform a joint test of HWD and sdMAF (i.e. `joint_test = TRUE`), `sdMAF` should be set to `NULL`. 
- `joint_test`: whether to jointly test HWD and sdMAF. Logical. If `joint_test = TRUE`, the function will perform a joint test of HWD and sdMAF. 

### Tutorial on how to test HWE for an NPR SNP

Suppose the genotype counts of an Xchr NPR SNP is $f_{AA} = 86$, $f_{AB} = 228$, $f_{BB} = 22$, $m_A = 222$, and $m_B = 94$. This examplary SNP is rs6655837 in the AFR super population of the 1000 Genome Project. 

```r
female_sample <- c(86, 228, 22)
male_sample <- c(222, 94)
```

### Joint test of HWE and no sdMAF

```r
RA_Xchr_github(gF = female_sample, gM = male_sample, snp_type='NPR', sdMAF = NULL, joint_test = TRUE)
[1] 1.06641e-15
```

### Test of HWE assuming no sdMAF

```r
RA_Xchr_github(gF = female_sample, gM = male_sample, snp_type='NPR', sdMAF = FALSE, joint_test = FALSE)
[1] 2.205595e-14
```

### Test of HWE assuming sdMAF

```r
RA_Xchr_github(gF = female_sample, gM = male_sample, snp_type='NPR', sdMAF = TRUE, joint_test = FALSE)
[1] 7.260256e-14
```


## Pseudoautosomal Regions (PARs)

In Xchr PAR, both males and females are diploid. Suppose the variant of interest is bi-allelic and A is the reference allele. Suppose we have a homogeneous population, female genotype counts are denoted as $f_{AA}$, $f_{AB}$ and $f_{BB}$; while male genotype counts are denoted as $m_{AA}$, $m_{AB}$ and $m_{BB}$. 

Use the `RA_Xchr` function to perform the test. Input arguments include
- `gF`: female genotype counts. A vector of length 3. `gF = c(fAA, fAB, fBB)`.
- `gM`: male genotype counts. A vector of length 3. `gM = c(mAA, mAB, mBB)`.
- `snp_type`: type of SNP. `snp_type = 'NPR' ` for Xchr NPR SNPs and `snp_type = 'PAR'` for Xchr PAR SNPs.
- `sdMAF`: sdMAF status. Logical. If `snp_type = TRUE`, the reference allele frequency is estimated separately in females and males; if `snp_type = FALSE`, the reference allele frequency is estimated using a pooled sample of females and males. 
- `joint_test`: whether to jointly test HWD and sdMAF. Logical. For a PAR SNP, `joint_test = FALSE` (no joint test of HWE and no sdMAF is available).

### Tutorial on how to test HWE for a PAR SNP

Suppose the genotype counts of an Xchr PAR SNP is $f_{AA} = 88$, $f_{AB} = 239$, $f_{BB} = 9$, $m_{AA} = 85$, $m_{AB} = 218$, and $m_{BB} = 13$. This examplary SNP is rs867436760 in the AFR super population of the 1000 Genome Project. 


```r
female_sample <- c(88, 239, 9)
male_sample <- c(85, 218, 13)
``` 

### Test of HWE assuming no sdMAF

```r
RA_Xchr_github(gF = female_sample, gM = male_sample, snp_type='PAR', sdMAF = FALSE, joint_test = FALSE)
[1] 1.032235e-34
```

### Test of HWE assuming sdMAF

```r
RA_Xchr_github(gF = female_sample, gM = male_sample, snp_type='PAR', sdMAF = TRUE, joint_test = FALSE)
[1] 1.273039e-33
```

