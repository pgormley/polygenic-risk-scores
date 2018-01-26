#!/bin/bash
#!/bin/bash
# Run like: qsub -o q.calc_PRS.NFBC66.out -l h_rt=06:00:00 -l h_vmem=4g ~/ugersub.sh sh calc_PRS_scores.sh NFBC66 ./imputation/NFBC/NFBC66 ./any_mig.nonFinns.gwama.out

# Name of study as prefix for output files e.g. swedish.QC4 || NFBC66 || GENMETS || LASERI
STUDY=$1
if [[ ! -d $STUDY ]]; then
	mkdir $STUDY
fi

# Root name of imputed data (in impute2/oxford format)
imputed_data=$2

# Select GWAS dataset (*must not include sample that will be scored*)
gwas_betas=$3

# Reference 1000Genomes data to be used for LD-clumping
LD_ref_data="../../1000g_plink_data/my.ALL_1000G_phase1integrated_v3_aug2012_macGT1.eur.bfile"
unrelateds="../../1000g_plink_data/keepUnrelated.txt"

# GWAS SNP filtering parameters to be used in PRS score
MAF=0.05	# GWAS filter SNPs by e.g. MAF > 5%
N_STUDIES=10	# GWAS filter SNPs by no. of studies successfully imputed e.g. N_STUDIES >= 10
I2=0.75		# GWAS filter SNPs by heterogeneity index e.g. i2 < 0.75

# Specify required column numbers in GWAS results file
snpID_column=3	# GWAS file column no. of SNP ID
allele_column=4	# GWAS file column no. of effect allele corresponding to beta
betas_column=7	# GWAS file column no. of estimated betas

# Filter GWAS betas on MAF, n_studies, and i2
awk -v maf=$MAF -v nstud=$N_STUDIES -v i2=$I2 'NR==1 || ($1!=-9 && $2!=-9 && $6>=maf && $6<=(1-maf) && $16<i2 && $17>=nstud){print $0}' $gwas_betas | sort -k1g,1 -k2g,2 > $gwas_betas.maf$MAF.isq$I2.nstud$N_STUDIES.clean

# Get list of imputed variants at INFO>0.6 and MAF>5% that intersect with list of filtered GWAS betas
for chr in {1..22}; do
	awk -v maf=0.05 -v info=0.6 'FNR==NR{if($4>=maf && $4<=(1-maf) && $5>=info){clean_snps[$2]=$2;next}else{next}}($3 in clean_snps || FNR==1){print $0}' ${imputed_data}*chr${chr}.*info $gwas_betas.maf$MAF.isq$I2.nstud$N_STUDIES.clean | sort -k1g,1 -k2g,2 > $STUDY/$STUDY.chr${chr}.intersection.gwama

	# Exclude duplicate SNP IDs in imputed dataset to avoid downstream errors when attempting to calculate PRS scores
	awk '{print $2}' ${imputed_data}*chr${chr}.*info | sort | uniq -d > $STUDY/$STUDY.chr${chr}.duplicate_snps
	grep -vwf $STUDY/$STUDY.chr${chr}.duplicate_snps $STUDY/$STUDY.chr${chr}.intersection.gwama > $STUDY/$STUDY.chr${chr}.intersection.gwama.no_dups
done

# LD-clump the intersecting variants list based on LD from 1000g EUR samples and p-values from GWAS betas file:
for chr in {1..22}; do
	plink --bfile $LD_ref_data --keep $unrelateds --maf 0.05 --clump $STUDY/$STUDY.chr${chr}.intersection.gwama.no_dups --clump-p1 0.05 --clump-p2 1.00 --clump-r2 0.1 --clump-kb 500 --clump-snp-field rs_number --clump-field p-value --out $STUDY/$STUDY.chr$chr
done

# Calculate PRS scores using list of LD-clumped SNPs:
for chr in {1..22}; do
	awk 'NR>1{print $3}' $STUDY/$STUDY.chr${chr}.clumped > $STUDY/$STUDY.chr${chr}.clumped.snplist
	#plink --gen ${imputed_data}*chr${chr}[_.]*gz --sample $imputed_data.sample --oxford-single-chr $chr --extract $STUDY/$STUDY.chr${chr}.clumped.snplist --score $STUDY/$STUDY.chr${chr}.intersection.gwama $snpID_column $allele_column $betas_column header --out $STUDY/$STUDY.chr$chr
	plink --gen ${imputed_data}*chr${chr}_*gz --sample $imputed_data.sample --oxford-single-chr $chr --extract $STUDY/$STUDY.chr${chr}.clumped.snplist --score $STUDY/$STUDY.chr${chr}.intersection.gwama $snpID_column $allele_column $betas_column header --out $STUDY/$STUDY.chr$chr
done

# Combine per-chromosome PRS scores by summing across the genome
cp $STUDY/$STUDY.chr1.profile $STUDY/$STUDY.profile.sum
for chr in {2..22}; do
	awk -v OFS="\t" 'FNR==NR{cnt[$1"-"$2]=$4; cnt2[$1"-"$2]=$5; score[$1"-"$2]=$6; next}($1"-"$2 in score){if(FNR==1){print $1,$2,$3,$4,$5,$6}else{print $1,$2,$3,$4+cnt[$1"-"$2],$5+cnt2[$1"-"$2],$6+score[$1"-"$2]}}' $STUDY/$STUDY.chr$chr.profile $STUDY/$STUDY.profile.sum > $STUDY/$STUDY.profile.sum.tmp
	mv $STUDY/$STUDY.profile.sum.tmp $STUDY/$STUDY.profile.sum
done
