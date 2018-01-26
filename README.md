# polygenic-risk-scores
----------------------
Author
----------------------
Padhraig Gormley

Massachusetts General Hospital | Broad Institute of MIT and Harvard.


----------------------
Synopsis
----------------------
Script for calculating polygenic risk scores in a genotyped sample using GWAS summary statistics from a disease/trait.  


----------------------
Dependencies
----------------------
PLINK (https://www.cog-genomics.org/plink2)


----------------------
Usage
----------------------
Open your bash or other shell. Then to run the script with the supplied toy data type:

sh ./calc-PRS-scores.sh myStudyName imputedSampleRootName GWASresultsFile
