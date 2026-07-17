#!/bin/bash

#----------------#
# vcfpipe v1.0.1 #
#----------------#

#--- Dependencies ---#
# vcftools
# bcftools
# htslib
# R

#-----------#
# Core code #
#-----------#

#--- module Help and definition of all flags ---#
usage='
vcfpipe version 1.0.1 Copyright (C) 2026 Lionel Di Santo

vcfpipe is a bash-based module of SeqForge designed for comprehensive filtering of VCF (Variant Call Format) files containing genetic variant (SNPs and indels) and invariant sites. It provides a unified framework to apply multiple layers of site and genotype filtering, as well as the detection of potential paralogous sites.

To run properly, vcfpipe needs the following dependencies:
- vcftools
- bcftools
- htslib
- R

** Note that dependencies are not provided with the module and must either be downloaded or loaded as modules by the user **

Before running the vcfpipe module, make sure that SeqForge and vcfpipe dependencies are added to your path (if installed locally) using the command lines below:
export PATH=/path_to_dependencies/Dependencies/:$PATH
export PATH=/path_to_SeqForge/SeqForge/:$PATH
!!! If dependencies are loaded as modules, you need not export any paths !!!

Reference (argument -m, --max-depth)
Li, H. (2014). Toward better understanding of artifacts in variant calling from high-coverage samples.
Bioinformatics, 30(20), 2843-2851. https://doi.org/10.1093/bioinformatics/btu356

Reference (argument -p, --paralogs)
McKinney, G.J., Waples, R.K., Seeb, L.W. and Seeb, J.E. (2017). Paralogs are revealed by proportion of heterozygotes and deviations in read ratios in 
genotyping-by-sequencing data from natural populations. Molecular Ecology Resources, 17, 656-669. https://doi.org/10.1111/1755-0998.12613

USAGE
SeqForge vcfpipe [options] <file.vcf(.gz)> <working_directory>
SeqForge vcfpipe version (to print the version of the software)

Options to be passed to the module:
    -h, --help                          Show this help page.
    -n, --missing-ind       <float>     This parameter filters the data according to per individual fraction of missing data (range between 0 and 1). All individuals within the provided VCF file
                                        with a fraction of missing data greater than this value will be removed. Default is 1 (keep all individuals).
    -q, --minQ              <integer>   Minimum site quality (QUAL). Include only sites with Quality value above this threshold (other sites are removed). Default value is 30.
    -c, --min-mac           <integer>   Minimum minor allele count (MAC). Include only sites with MAC greater than or equal to this value (other sites are removed). Default value is 3.
    -g, --minGQ             <integer>   Minimum genotype quality (GQ). Exclude all genotypes with quality below this threshold (GQ FORMAT tag must be specified for all sites).
                                        Note that excluded genotypes are treated as missing (and not discarded). Default value is 20.
    -d, --minDP             <integer>   Minimum genotype depth (MinDP). Include only genotypes with depth greater than or equal to this value (DP FORMAT tag must be specified for all sites).
                                        Note that excluded genotypes are treated as missing (and not discarded). Default value is 5.
    -a, --maf               <float>     Minor allele frequency (MAF). Include only sites with a MAF greater than or equal to this value (other sites are removed). Allele frequency is defined 
                                        as the number of times an allele appears over all individuals at that site, divided by the total number of non-missing alleles at that site. Default value is 0.01.
    -i, --GT-call           <float>     Genotype call rate across all individuals. Exclude sites on the basis of the proportion of missing data (range between 0 and 1, where 0 allows sites
                                        that are completely missing, and 1 indicates no missing data allowed). Sites not complying with this argument are removed. Note that by default the module 
                                        will exlude any sites with a proportion of missing data greater than 0.5 before implementing this filter. Consequently, the value provided here should be 
                                        greater than 0.5 to have any incidence on retained sites. Default value is 0.95.
    -l, --linkage                       If specified (-l, --linkage), linkage desequilibrium filtering will be implemented and -r (--r2) and -s (--wd-size) must be provided (see below).
    -r, --r2                <float>     R2 for linkage desequilibrium filtering. This represents the maximum correlation squared allow between two genetic variants. Any genetic variants with 
                                        a correlation squared above this threshold (within the window specified by -w, --wd-size (see below) will be removed.
    -s, --wd-size           <integer>   Window size used during linkage desequilibrium filtering (in bp). For more details, see above (-r, --r2).
    -z, --min-depth         <integer>   Minimum average read depth (over all individuals). Include only sites with a mean read depth over all individual greater than or equal to this value. Default valus is 5.
    -m, --max-depth         <string>    This argument takes one of the following values: Li2014, mean:x, none, or a numerical value. If 'Li2014' is provided (-m, --max-depth Li2014), then filtering of sites 
                        or  <integer>   based on maximum read depth following recommendations by Li 2014 is performed (see above for full reference); sites with mean depth (over all individuals) > d+4*sqrt(d),
                                        where d is the average read depth across variants, are removed. If 'mean:x' is provided (e.g., -m, --max-depth mean:2), where x must be a numerical value, then sites with
                                        mean read depth (over all individuals) > x times the mean read depth across variants are removed. If 'none' is provided (-m, --max-depth none), no filtering based on site
                                        mean depth will be performed. If a numerical value is provided (e.g., -m, --max-depth 100), then sites with mean read depth (over all individuals) > than the provided value
                                        (e.g., 100) will be removed. Default is 'none'.
    -p, --paralogs          <string>    This argument takes one of the following values: none or x:y. If 'x:y' is provided, genetic variants likely stemming from the mapping of paralogous sequences will be identified
                                        and removed from the VCF file. Identification of paralogs necessitates two parameters, namely x and y. x fixes the upper threshold for the proportion of heterozygotes 
                                        [H, ranges from >0 to 1]; any sites with heterozygosity above this threshold will be discarded. y fixes the upper limit for the deviation of ratios from the expected 50:50 
                                        [D, range from -inf to inf], calcuated as a z-score. In other words, this parameter defines how many standard deviations the observed A:B allele ratio at a given heterozygous
                                        site is allowed to deviate from the expected 50:50 balance. This parameter thus filters sites with strong allelic imbalance. Note that when using this parameter, it is best to
                                        use only bi-allelic genetic markers. This option supports multi-allelic genetic markers, but sites with more than two alleles are ignored during heterozygosity and allelic imbalance
                                        filtering. For more details on this approach and how to select H and D thresholds, please refer to the original paper by McKinney et al. 2017 (see above for the full reference). 
                                        If 'none' is provided, no filtering of paralogs will be performed. To use this option, R package vcfR, ggplot2, dplyr, and stringr must be installed. Double-check that the packages
                                        are installed for the R version loaded or exported. If an R repository containing vcfR, ggplot2, dplyr, and stringr packages is missing, the module will return an error.
                                        Default is 'none'.
    -t, --threads           <integer>   Set the number of threads to use for analysis. Default is 1. 
    -b, --no-multiallelic               If specified (-b, --no-multiallelic), this argument triggers the exclusion of multiallelic sites, so only biallelic single-nucleotide polymorphisms and 
                                        indels are kept.
    -e, --no-indels                     If specified (-e, --no-indels), this argument triggers the exclusion of indels, so that exclusively single-nucleotide polymorphisms are kept (either biallelic 
                                        or multiallelic depending on whether -b, --no-multiallelic is specified).
'
if [ $# -eq 0 ]
then
    echo ""
    echo "$usage"
    echo ""
    exit 0
fi

PARSED_OPTIONS=$(getopt -o hn:q:c:g:d:a:i:lr:s:z:m:p:t:be -l help,missing-ind:,minQ:,min-mac:,minGQ:,minDP:,maf:,GT-call:,linkage,r2:,wd-size:,min-depth:,max-depth:,paralogs:,threads:,no-multiallelic,no-indels -- "$@")
if [[ $? -ne 0 ]]; then echo -e "\nError with parsing arguments\n"; exit 1; fi
eval set -- "$PARSED_OPTIONS"

while true; do
    case "$1" in
        -h|--help) 
            echo "$usage"
            exit
            ;;
        -n|--missing-ind)
            indmiss="$2"
            shift 2
            ;;
        -q|--minQ)
            minq="$2"
            shift 2
            ;;
        -c|--min-mac)
            mac="$2"
            shift 2
            ;;
        -g|--minGQ)
            mingq="$2"
            shift 2
            ;;
        -d|--minDP)
            mindp="$2"
            shift 2
            ;;
        -a|--maf)
            maf="$2"
            shift 2
            ;;
        -i|--GT-call)
            gcall="$2"
            shift 2
            ;;
        -l|--linkage)
            ld=yes
            shift
            ;;
        -r|--r2)
            rsq="$2"
            shift 2
            ;;
        -s|--wd-size)
            win="$2"
            shift 2
            ;;
        -z|--min-depth)
            minmeandp="$2"
            shift 2
            ;;
        -m|--max-depth)
            meandp="$2"
            shift 2
            ;;
        -p|--paralogs)
            paralogs="$2"
            shift 2
            ;; 
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -b|--no-multiallelic)
            biallelic=yes
            shift
            ;;
        -e|--no-indels)
            indels=yes
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Unexpected option: $1"
            exit 1
            ;;
    esac
done

if [ $# -eq 1 ] && [ $1 = version ]
    then
        echo ""
        echo "vcfpipe version 1.0.1 Copyright (C) 2026 Lionel Di Santo"
        echo ""
        exit 0
fi
if [ $# -lt 2 ]
    then
        echo ""
        echo "The VCF file and/or the working directory is/are missing."
        echo ""
        exit 1
fi
if [ $# -gt 2 ]
    then
        echo ""
        echo "Too many arguments provided! Two arguments only must be provided following [options]. The VCF file and the working directory."
        echo ""
        exit 1
fi
if [ ! -z $ld ] && [ -z $rsq ]
    then
        echo ""
        echo "When using -l (--linkage), an R2 and window size value needs to be provided using -r (--r2) and -s (--wd-size) respectively (see SeqForge vcfpipe -h for help)."
        echo ""
        exit 1
fi
if [ ! -z $ld ] && [ -z $win ]
    then
        echo ""
        echo "When using -l (--linkage), an R2 and window size value needs to be provided using -r (--r2) and -s (--wd-size) respectively (see SeqForge vcfpipe -h for help)."
        echo ""
        exit 1
fi

#--- Setting the working directory ---#
cd $2

#--- Welcome message ---#
echo ""
echo "#----------------------------------------------------------#"
echo "# vcfpipe version 1.0.1 Copyright (C) 2026 Lionel Di Santo #"
echo "#----------------------------------------------------------#"
program=("R" "vcftools" "bcftools" "htsfile")
for i in ${program[@]}
    do
        if command -v $i > /dev/null 2>&1
            then
                echo "$(basename $i | sed 's/file/lib/') is installed."
            else 
                echo "$(basename $i | sed 's/file/lib/') is not installed or cannot be found. Please make sure to install $(basename $i | sed 's/file/lib/') or export it to your path."
                ext=true
        fi
done
if [ $ext ]; then exit 1; fi

echo ""
echo "Parameters defined for vcfpipe run"
echo "----------------------------------"

#--- Print run parameters ---#
if [ -z $threads ]
    then
        echo "Number of threads for compression: 1"
        threads=1
    else
        echo "Number of threads for compression: $threads"
fi
if [ ! -z $ld ]
    then
        echo "Option -l, --linkage activated: Genetic variants will be filtered for linkage desequilibrium."
fi
if [ ! -z $ld ]
    then
        echo "R2 for linkage desequilibrium filtering: $rsq"
fi
if [ ! -z $ld ]
    then
        echo "Window size used during linkage desequilibrium filtering: $win bp"
fi
if [ ! -z $biallelic ]
    then
        echo "Option -b, --no-multiallelic activated: Only biallelic single-nucleotide polymorphisms and indels will be kept in the vcf file (multiallelic SNPs are removed)."
fi
if [ ! -z $indels ]
    then
        echo "Option -e, --no-indels activated: Only single-nucleotide polymorphisms will be kept in the vcf file (indels are removed). SNPs might be either biallelic or multiallelic depending on whether -b (--no-multiallelic) is specified."
fi
if [ -z $indmiss ]
    then
        echo "All individuals kept (no filtering based on per individual missingness)."
        indmiss='1'
    else 
        indmiss="$indmiss"
        echo "Maximum fraction of missing data allowed to keep an individual: $indmiss"
fi
if [ -z $minmeandp ]
    then
        echo "Sites with mean read depth (over all individuals) < 5 will be removed."
        minmeandp=5
    else 
        echo "Sites with mean read depth (over all individuals) < $minmeandp will be removed."
fi
if [ -z $meandp ] 
    then
        meandp=none
    else
        if [[ $meandp != none ]]
            then
                if [[ $meandp = Li2014 ]]; then echo "Sites will be filtered based on mean read depth (over all individuals) following recommendations by Li 2014."; fi
                if [[ $meandp =~ ^mean:[0-9]+$ ]]; then echo "Sites with mean read depth (over all individuals) > $(echo $meandp | cut -f2 -d:) times the mean read depth across variants will be removed."; fi
                if [[ $meandp =~ ^[0-9]+$ ]]; then echo "Sites with mean read depth (over all individuals) > $meandp will be removed."; fi
                if [[ ! $meandp = Li2014 ]] && [[ ! $meandp = none ]] && [[ ! $meandp =~ ^mean:[0-9]+$ ]] && [[ ! $meandp =~ ^[0-9]+$  ]]; then echo ""; echo "String provided for -m, --max-depth is invalid (see SeqForge vcfpipe -h for help)."; echo ""; exit 1; fi
        fi
fi
if [ -z $minq ]
    then
        echo "Minimum site quality (QUAL): 30"
        minq=30
    else
        echo "Minimum site quality (QUAL): $minq"
fi
if [ -z $mac ]
    then
        echo "Minimum allele count (MAC): 3"
        mac=3
    else
        echo "Minimum allele count (MAC): $mac"
fi
if [ -z $mingq ]
    then
        echo "Minimum genotype quality (GQ): 20"
        mingq=20
    else
        echo "Minimum genotype quality (GQ): $mingq"
fi
if [ -z $mindp ]
    then
        echo "Minimum genotype depth (MinDP): 5"
        mindp=5
    else
        echo "Minimum genotype depth (MinDP): $mindp"
fi
if [ -z $maf ]
    then
        echo "Minor allele frequency (MAF): 0.01"
        maf=0.01
    else
        echo "Minor allele frequency (MAF): $maf"
fi
if [ -z $gcall ]
    then
        echo "Genotype call rate across all individuals: 0.95"
        gcall=0.95
    else
        echo "Genotype call rate across all individuals: $gcall"
fi
if [ -z $paralogs ] 
    then
        paralogs=none
    else
        if [[ $paralogs != none ]]
            then
                if [[ ! $paralogs = none ]] && [[ ! $paralogs =~ ^([0-9]\.[0-9]+:)?(1:)?[0-9]+$ ]]; then echo ""; echo "String provided for -p, --paralogs is invalid (see SeqForge vcfpipe -h for help)."; echo ""; exit 1; fi
                het=$(echo $paralogs | cut -f1 -d:); std=$(echo $paralogs | cut -f2 -d:)
                echo "Once all other filters are applied, paralogous sites will be identified and discarded based on H > $het and |D| > $std."
        fi
fi
echo ""
echo "VCF file provided: $1"
echo "Working directory provided: $(pwd)"
echo ""

#--- Filtering of genetic variants ---#
echo "Filtering of genetic variants"
echo "-----------------------------"

if [[ $1 = *.gz ]]
    then
        echo "VCF file provided is compressed - Uncompressing ${1}..."
        bgzip -@ $threads -d $1
        inputVCF=${1%.gz}
        echo ""
    else
        inputVCF=$1
fi

if [ -e vcftools.log ]; then rm vcftools.log; fi
if [ -e bcftools.log ]; then rm bcftools.log; fi
check=$(bcftools view -i 'ALT="."' $inputVCF | grep -v '^#' | wc -l)
if [[ $check -eq 0 ]]
    then
        echo "Marking genotypes with < $mindp reads and a quality < $mingq as missing"
        echo "  & removing sites with more than 50% missing genotypes..."
        vcftools --vcf $inputVCF --minGQ $mingq --minDP $mindp --max-missing 0.5 --recode --out SNP_temp_missing >> vcftools.log 2>&1
        echo ""

        if [ $indmiss != '1' ]
            then
                echo "Removing individuals with fraction of missing data > $indmiss..."
                vcftools --vcf SNP_temp_missing.recode.vcf --missing-indv --out ind >> vcftools.log 2>&1
                echo 'args = commandArgs(trailingOnly = TRUE)' > indFilt.R
                echo 'miss <- read.table("ind.imiss", h = T)' >> indFilt.R
                echo 'nrow(miss)' >> indFilt.R
                echo 'range(miss$F_MISS)' >> indFilt.R
                echo 'pdf(file = "Fractions_missing_individuals.pdf", width = 9, height = 6)' >> indFilt.R
                echo 'hist(miss$F_MISS, main = "", xlab = "Fraction of missing data (per individual)")' >> indFilt.R
                echo 'dev.off()' >> indFilt.R
                echo 'max.missing <- as.numeric(args[1])' >> indFilt.R
                echo 'max.missing' >> indFilt.R
                echo 'keep <- miss[which(miss$F_MISS <= max.missing),]' >> indFilt.R
                echo 'write.table(keep$INDV, file = "inds.keep", quote = F, row.names = F, col.names = F, sep = "\n")' >> indFilt.R
                echo 'sessionInfo()' >> indFilt.R
                echo 'sessioninfo::package_info()' >> indFilt.R
                R CMD BATCH --vanilla "--args $indmiss" indFilt.R
                vcftools --vcf SNP_temp_missing.recode.vcf --keep inds.keep --recode --out indFiltSNPs >> vcftools.log 2>&1
                echo ""
        fi

        if [ $indmiss != '1' ]
            then
                vcF=indFiltSNPs.recode.vcf
            else
                vcF=SNP_temp_missing.recode.vcf
        fi

        echo "Removing sites from VCF file with:" 
        echo "  - minimum allele count < $mac"
        echo "  - minimum site quality <= $minq"
        echo "  - minimum read depth over all individuals < $minmeandp"
        echo "  - minor allele frequency < $maf"
        echo "  - Genotype call rate < $gcall" 
        echo "  Using $vcF for filtering..."
        vcftools --vcf $vcF --minQ $minq --mac $mac --min-meanDP $minmeandp --max-missing $gcall --maf $maf --recode --out SNP_temp >> vcftools.log 2>&1
        echo ""

        if [[ $meandp != none ]]
            then
                if [[ $meandp =~ ^mean:[0-9]+$ ]]; then xmean=$(echo $meandp | cut -f2 -d:); echo "Removing variants with mean depth > $xmean times mean read depth across variants from VCF file..."; fi
                if [[ $meandp = Li2014 ]]; then echo "Removing variants with mean depth > d+4*sqrt(d) from VCF file..."; fi
                if [[ ! $meandp =~ ^[0-9]+$ ]]
                    then
                        echo "  Estimating mean depth..."
                        vcftools --vcf $inputVCF --site-mean-depth --out max_meanDP >> vcftools.log 2>&1
                        echo 'df <- read.table("max_meanDP.ldepth.mean", h=T)' > max_meanDP.R
                        echo 'mean <- mean(df[,3]); mean' >> max_meanDP.R
                        if [[ $meandp = Li2014 ]]; then echo 'res <- ceiling(mean+(4*sqrt(mean))); res' >> max_meanDP.R; echo 'write(res, file = "max_meanDP_estimate.txt")' >> max_meanDP.R; fi
                        if [[ $meandp =~ ^mean:[0-9]+$ ]]; then echo "xmean <- round($xmean*mean,0); xmean" >> max_meanDP.R; echo 'write(xmean, file = "max_meanDP_estimate.txt")' >> max_meanDP.R; fi
                        echo 'sessionInfo()' >> max_meanDP.R
                        echo 'sessioninfo::package_info()' >> max_meanDP.R
                        R CMD BATCH --vanilla max_meanDP.R
                        max_mean_DP_est=$(cat max_meanDP_estimate.txt)
                    else
                        max_mean_DP_est=$meandp
                fi
                echo "  Removal of variants with mean depth > $max_mean_DP_est"
                vcftools --vcf SNP_temp.recode.vcf --max-meanDP $max_mean_DP_est --recode --out SNP_temp2 >> vcftools.log 2>&1
                echo ""
        fi

        if [ ! -z $biallelic ]
            then
                echo "Removing multiallelic sites from VCF file..."
                if [ -e SNP_temp2.recode.vcf ]; then file=SNP_temp2.recode.vcf; else file=SNP_temp.recode.vcf; fi
                echo "  Using $file for filtering..."
                vcftools --vcf $file --max-alleles 2 --recode --out SNP_temp3 >> vcftools.log 2>&1
                echo ""
        fi

        if [ ! -z $indels ]
            then
                echo "Removing indels from VCF file..."
                if [ -e SNP_temp2.recode.vcf ]; then file=SNP_temp2.recode.vcf; fi
                if [ -e SNP_temp3.recode.vcf ]; then file=SNP_temp3.recode.vcf; fi
                list=("SNP_temp2.recode.vcf" "SNP_temp3.recode.vcf")
                for j in ${list[@]}
                    do
                        if [ -e $j ]; then any_file=true; break; fi
                done
                if [ ! $any_file ]; then file=SNP_temp.recode.vcf; fi
                echo "  Using $file for filtering..."
                vcftools --vcf $file --remove-indels --recode --out SNP_temp4 >> vcftools.log 2>&1
                echo ""
        fi

        if [ ! -z $ld ]
            then
                echo "Filtering sites based on linkage desequilibrium (R2 = $rsq) within a window of $win bp..."
                if [ -e SNP_temp2.recode.vcf ]; then file=SNP_temp2.recode.vcf; fi
                if [ -e SNP_temp3.recode.vcf ]; then file=SNP_temp3.recode.vcf; fi
                if [ -e SNP_temp4.recode.vcf ]; then file=SNP_temp4.recode.vcf; fi
                list=("SNP_temp2.recode.vcf" "SNP_temp3.recode.vcf" "SNP_temp4.recode.vcf")
                for j in ${list[@]}
                    do
                        if [ -e $j ]; then any_file=true; break; fi
                done
                if [ ! $any_file ]; then file=SNP_temp.recode.vcf; fi
                echo "  Using $file for filtering ..."
                bgzip -@ $threads $file
                bcftools index --threads $threads ${file}.gz
                bcftools +prune -m $rsq -w ${win}bp -Ov -o SNP_temp5.vcf ${file}.gz >> bcftools.log 2>&1
                echo ""
        fi

        if [[ $paralogs != none ]]
            then
                het=$(echo $paralogs | cut -f1 -d:)
                std=$(echo $paralogs | cut -f2 -d:)
                echo 'Identification and removal of paralogous sites...'
                echo "  - H limit set to $het"
                echo "  - |D| limit set to $std"
                if [[ -e HDplots.R ]]; then echo "Previous run with --paralogs detected and discarded."; rm HDplots.R*; fi
                bash create_HDplots.sh
                if [ -e SNP_temp2.recode.vcf ]; then file=SNP_temp2.recode.vcf; fi
                if [ -e SNP_temp3.recode.vcf ]; then file=SNP_temp3.recode.vcf; fi
                if [ -e SNP_temp4.recode.vcf ]; then file=SNP_temp4.recode.vcf; fi
                if [ -e SNP_temp5.vcf ]; then file=SNP_temp5.vcf; fi
                list=("SNP_temp2.recode.vcf" "SNP_temp3.recode.vcf" "SNP_temp4.recode.vcf" "SNP_temp5.vcf")
                for j in ${list[@]}
                    do
                        if [ -e $j ]; then any_file=true; break; fi
                done
                if [ ! $any_file ]; then file=SNP_temp.recode.vcf; fi
                echo "  Using $file for the identification and removal of paralogous sites ..."
                R CMD BATCH --vanilla "--args $file $het $std" HDplots.R
                if [[ ! -e keep_paralogs.txt ]]; then echo ""; echo "Something went wrong while identifying paralogous site. It may be that not all require R packages are installed. Check HDplots.Rout for details"; echo ""; exit 1; fi
                if [[ $(cat keep_paralogs.txt | wc -l) -gt 0 ]]
                    then
                        bgzip -@ $threads $file
                        bcftools sort -o ${file%.vcf}.sorted.vcf.gz -T . ${file}.gz > /dev/null 2>&1
                        bcftools index --threads $threads ${file%.vcf}.sorted.vcf.gz
                        bcftools view -R keep_paralogs.txt -Oz -o ${file%.vcf}.FiltD.vcf.gz ${file%.vcf}.sorted.vcf.gz >> bcftools.log 2>&1
                        mv ${file%.vcf}.FiltD.vcf.gz FinalSNPs.vcf.gz
                        if [ -e indFiltSNPs.recode.vcf ]; then rm indFiltSNPs.recode.vcf; fi
                        rm SNP_temp*
                    else
                        bgzip -@ $threads $file
                        mv $file FinalSNPs.vcf.gz
                        if [ -e indFiltSNPs.recode.vcf ]; then rm indFiltSNPs.recode.vcf; fi
                        rm SNP_temp*
                fi
        fi

        if [[ ! -e FinalSNPs.vcf.gz ]]
            then
                if [ -e indFiltSNPs.recode.vcf ]; then rm indFiltSNPs.recode.vcf; fi
                if [ -e SNP_temp5.vcf ]; then mv SNP_temp5.vcf FinalSNPs.vcf; rm SNP_temp*; fi
                if [ -e SNP_temp4.recode.vcf ]; then mv SNP_temp4.recode.vcf FinalSNPs.vcf; rm SNP_temp*; fi
                if [ -e SNP_temp3.recode.vcf ]; then mv SNP_temp3.recode.vcf FinalSNPs.vcf; rm SNP_temp*; fi
                if [ -e SNP_temp2.recode.vcf ]; then mv SNP_temp2.recode.vcf FinalSNPs.vcf; rm SNP_temp*; fi
                if [ -e SNP_temp.recode.vcf ]; then mv SNP_temp.recode.vcf FinalSNPs.vcf; rm SNP_temp*; fi
                bgzip -@ $threads FinalSNPs.vcf
        fi
        bgzip -@ $threads $inputVCF

    else
        echo "NOTE: Invariant sites detected in the VCF file - $inputVCF (number = $check)!"
        echo "Filters on minor alleles (--min-mac, --maf, and --paralogs) will be applied to variant sites only".
        echo "Multiple output VCF will be produced:"
        echo "  - FinalSNPs_invariant.vcf.gz: VCF file containing only invariant sites filtered for all options but --min-mac $mac, --maf $maf, and --paralogs $paralogs."
        echo "  - FinalSNPs_variant.vcf.gz: VCF file containing only variant sites filtered for all options, including --min-mac $mac, --maf $maf, and --paralogs $paralogs."
        echo "  - FinalSNPs_variant_invariant.vcf.gz: A file concatenating sites present in FinalSNPs_invariant.vcf.gz and FinalSNPs_variant.vcf.gz."
        echo ""

        echo "Marking genotypes with < $mindp reads and a quality < $mingq as missing" 
        echo "  & removing sites with more than 50% missing genotypes..."
        vcftools --vcf $inputVCF --minGQ $mingq --minDP $mindp --max-missing 0.5 --recode --out SNP_temp_missing >> vcftools.log 2>&1
        echo ""

        if [ $indmiss != '1' ]
            then
                echo "Removing individuals with fraction of missing data > $indmiss..."
                vcftools --vcf SNP_temp_missing.recode.vcf --missing-indv --out ind >> vcftools.log 2>&1
                echo 'args = commandArgs(trailingOnly = TRUE)' > indFilt.R
                echo 'miss <- read.table("ind.imiss", h = T)' >> indFilt.R
                echo 'nrow(miss)' >> indFilt.R
                echo 'range(miss$F_MISS)' >> indFilt.R
                echo 'pdf(file = "Fractions_missing_individuals.pdf", width = 9, height = 6)' >> indFilt.R
                echo 'hist(miss$F_MISS, main = "", xlab = "Fraction of missing data (per individual)")' >> indFilt.R
                echo 'dev.off()' >> indFilt.R
                echo 'max.missing <- as.numeric(args[1])' >> indFilt.R
                echo 'max.missing' >> indFilt.R
                echo 'keep <- miss[which(miss$F_MISS <= max.missing),]' >> indFilt.R
                echo 'write.table(keep$INDV, file = "inds.keep", quote = F, row.names = F, col.names = F, sep = "\n")' >> indFilt.R
                echo 'sessionInfo()' >> indFilt.R
                echo 'sessioninfo::package_info()' >> indFilt.R
                R CMD BATCH --vanilla "--args $indmiss" indFilt.R
                vcftools --vcf SNP_temp_missing.recode.vcf --keep inds.keep --recode --out indFiltSNPs >> vcftools.log 2>&1
                echo ""
        fi

        if [ $indmiss != '1' ]
            then
                vcF=indFiltSNPs.recode.vcf
            else
                vcF=SNP_temp_missing.recode.vcf
        fi

        echo "Removing sites from VCF file with:" 
        echo "  - minimum site quality <= $minq"
        echo "  - minimum read depth over all individuals < $minmeandp"
        echo "  - Genotype call rate < $gcall" 
        echo "  Using $vcF for filtering..."
        vcftools --vcf $vcF --minQ $minq --min-meanDP $minmeandp --max-missing $gcall --recode --out SNP_temp >> vcftools.log 2>&1
        echo ""

        if [[ $meandp != none ]]
            then
                if [[ $meandp =~ ^mean:[0-9]+$ ]]; then xmean=$(echo $meandp | cut -f2 -d:); echo "Removing variants with mean depth > $xmean times mean read depth across variants from VCF file..."; fi
                if [[ $meandp = Li2014 ]]; then echo "Removing variants with mean depth > d+4*sqrt(d) from VCF file..."; fi
                if [[ ! $meandp =~ ^[0-9]+$ ]]
                    then
                        echo "  Estimating mean depth..."
                        vcftools --vcf $inputVCF --site-mean-depth --out max_meanDP >> vcftools.log 2>&1
                        echo 'df <- read.table("max_meanDP.ldepth.mean", h=T)' > max_meanDP.R
                        echo 'mean <- mean(df[,3]); mean' >> max_meanDP.R
                        if [[ $meandp = Li2014 ]]; then echo 'res <- ceiling(mean+(4*sqrt(mean))); res' >> max_meanDP.R; echo 'write(res, file = "max_meanDP_estimate.txt")' >> max_meanDP.R; fi
                        if [[ $meandp =~ ^mean:[0-9]+$ ]]; then echo "xmean <- round($xmean*mean,0); xmean" >> max_meanDP.R; echo 'write(xmean, file = "max_meanDP_estimate.txt")' >> max_meanDP.R; fi
                        echo 'sessionInfo()' >> max_meanDP.R
                        echo 'sessioninfo::package_info()' >> max_meanDP.R
                        R CMD BATCH --vanilla max_meanDP.R
                        max_mean_DP_est=$(cat max_meanDP_estimate.txt)
                    else
                        max_mean_DP_est=$meandp
                fi
                echo "  Removal of variants with mean depth > $max_mean_DP_est"
                vcftools --vcf SNP_temp.recode.vcf --max-meanDP $max_mean_DP_est --recode --out SNP_temp2 >> vcftools.log 2>&1
                echo ""
        fi

        if [ ! -z $biallelic ]
            then
                echo "Removing multiallelic sites from VCF file..."
                if [ -e SNP_temp2.recode.vcf ]; then file=SNP_temp2.recode.vcf; else file=SNP_temp.recode.vcf; fi
                echo "  Using $file for filtering..."
                vcftools --vcf $file --max-alleles 2 --recode --out SNP_temp3 >> vcftools.log 2>&1
                echo ""
        fi

        if [ ! -z $indels ]
            then
                echo "Removing indels from VCF file..."
                if [ -e SNP_temp2.recode.vcf ]; then file=SNP_temp2.recode.vcf; fi
                if [ -e SNP_temp3.recode.vcf ]; then file=SNP_temp3.recode.vcf; fi
                list=("SNP_temp2.recode.vcf" "SNP_temp3.recode.vcf")
                for j in ${list[@]}
                    do
                        if [ -e $j ]; then any_file=true; break; fi
                done
                if [ ! $any_file ]; then file=SNP_temp.recode.vcf; fi
                echo "  Using $file for filtering..."
                vcftools --vcf $file --remove-indels --recode --out SNP_temp4 >> vcftools.log 2>&1
                echo ""
        fi

        if [ ! -z $ld ]
            then
                echo "Filtering sites based on linkage desequilibrium (R2 = $rsq) within a window of $win bp..."
                if [ -e SNP_temp2.recode.vcf ]; then file=SNP_temp2.recode.vcf; fi
                if [ -e SNP_temp3.recode.vcf ]; then file=SNP_temp3.recode.vcf; fi
                if [ -e SNP_temp4.recode.vcf ]; then file=SNP_temp4.recode.vcf; fi
                list=("SNP_temp2.recode.vcf" "SNP_temp3.recode.vcf" "SNP_temp4.recode.vcf")
                for j in ${list[@]}
                    do
                        if [ -e $j ]; then any_file=true; break; fi
                done
                if [ ! $any_file ]; then file=SNP_temp.recode.vcf; fi
                echo "  Using $file for filtering ..."
                bgzip -@ $threads $file
                bcftools index --threads $threads ${file}.gz
                bcftools +prune -m $rsq -w ${win}bp -Ov -o SNP_temp5.vcf ${file}.gz >> bcftools.log 2>&1
                echo ""
        fi

        echo "Splitting input VCF into two files (one containing only invariant sites and another containing only variant sites)..."
        if [ -e SNP_temp2.recode.vcf ]; then file=SNP_temp2.recode.vcf; fi
        if [ -e SNP_temp3.recode.vcf ]; then file=SNP_temp3.recode.vcf; fi
        if [ -e SNP_temp4.recode.vcf ]; then file=SNP_temp4.recode.vcf; fi
        if [ -e SNP_temp5.vcf ]; then file=SNP_temp5.vcf; fi
        list=("SNP_temp2.recode.vcf" "SNP_temp3.recode.vcf" "SNP_temp4.recode.vcf" "SNP_temp5.vcf")
        for j in ${list[@]}
            do
                if [ -e $j ]; then any_file=true; break; fi
        done
        if [ ! $any_file ]; then file=SNP_temp.recode.vcf; fi        
        echo "  Using $file for splitting..."
        vcftools --vcf $file --max-maf 0 --recode --out invariant >> vcftools.log 2>&1
        vcftools --vcf $file --mac 1 --recode --out variant >> vcftools.log 2>&1
        echo ""

        echo "Removing (variant) sites from VCF file with:" 
        echo "  - minimum allele count < $mac"
        echo "  - minor allele frequency < $maf"
        vcftools --vcf variant.recode.vcf --mac $mac --maf $maf --recode --out variant_mac_maf >> vcftools.log 2>&1
        echo ""

        if [[ $paralogs != none ]]
            then
                het=$(echo $paralogs | cut -f1 -d:)
                std=$(echo $paralogs | cut -f2 -d:)
                echo 'Identification and removal of paralogous (variant) sites...'
                echo "  - H limit set to $het"
                echo "  - |D| limit set to $std"
                echo ""
                if [[ -e HDplots.R ]]; then echo "Previous run with --paralogs detected and discarded."; rm HDplots.R*; fi
                bash create_HDplots.sh
                file=variant_mac_maf.recode.vcf
                R CMD BATCH --vanilla "--args $file $het $std" HDplots.R
                if [[ ! -e keep_paralogs.txt ]]; then echo ""; echo "Something went wrong while identifying paralogous site. It may be that not all require R packages are installed. Check HDplot.Rout for details"; echo ""; exit 1; fi
                echo 'Concatenation of variant and invariant VCF files...'
                if [[ $(cat keep_paralogs.txt | wc -l) -gt 0 ]]
                    then
                        bgzip -@ $threads variant_mac_maf.recode.vcf
                        bgzip -@ $threads invariant.recode.vcf
                        bcftools sort -o variant_mac_maf_sorted.vcf.gz -T . variant_mac_maf.recode.vcf.gz > /dev/null 2>&1
                        bcftools sort -o invariant_sorted.vcf.gz -T . invariant.recode.vcf.gz > /dev/null 2>&1
                        bcftools index --threads $threads variant_mac_maf_sorted.vcf.gz
                        bcftools index --threads $threads invariant_sorted.vcf.gz
                        bcftools view -R keep_paralogs.txt -Oz -o variant_mac_maf_FiltD.vcf.gz variant_mac_maf_sorted.vcf.gz >> bcftools.log 2>&1
                        bcftools sort -o variant_mac_maf_FiltD_sorted.vcf.gz -T . variant_mac_maf_FiltD.vcf.gz > /dev/null 2>&1
                        bcftools index --threads $threads variant_mac_maf_FiltD_sorted.vcf.gz
                        bcftools concat --threads $threads --allow-overlaps -Oz -o FinalSNPs_variant_invariant.vcf.gz invariant_sorted.vcf.gz variant_mac_maf_FiltD_sorted.vcf.gz >> bcftools.log 2>&1
                        mv variant_mac_maf_FiltD_sorted.vcf.gz FinalSNPs_variant.vcf.gz
                        mv invariant_sorted.vcf.gz FinalSNPs_invariant.vcf.gz
                        rm variant* invariant*
                    else
                        bgzip -@ $threads variant_mac_maf.recode.vcf
                        bgzip -@ $threads invariant.recode.vcf
                        bcftools sort -o invariant_sorted.vcf.gz -T . invariant.recode.vcf.gz > /dev/null 2>&1
                        bcftools sort -o variant_mac_maf_sorted.vcf.gz -T . variant_mac_maf.recode.vcf.gz > /dev/null 2>&1
                        bcftools index --threads $threads invariant_sorted.vcf.gz
                        bcftools index --threads $threads variant_mac_maf_sorted.vcf.gz
                        bcftools concat --threads $threads --allow-overlaps -Oz -o FinalSNPs_variant_invariant.vcf.gz invariant_sorted.vcf.gz variant_mac_maf_sorted.vcf.gz >> bcftools.log 2>&1
                        mv variant_mac_maf.recode.vcf.gz FinalSNPs_variant.vcf.gz
                        mv invariant.recode.vcf.gz FinalSNPs_invariant.vcf.gz
                        rm variant* invariant*
                fi
            else
                echo 'Concatenation of variant and invariant VCF files...'
                bgzip -@ $threads variant_mac_maf.recode.vcf
                bgzip -@ $threads invariant.recode.vcf
                bcftools sort -o invariant_sorted.vcf.gz -T . invariant.recode.vcf.gz > /dev/null 2>&1
                bcftools sort -o variant_mac_maf_sorted.vcf.gz -T . variant_mac_maf.recode.vcf.gz > /dev/null 2>&1
                bcftools index --threads $threads invariant_sorted.vcf.gz
                bcftools index --threads $threads variant_mac_maf_sorted.vcf.gz
                bcftools concat --threads $threads --allow-overlaps -Oz -o FinalSNPs_variant_invariant.vcf.gz invariant_sorted.vcf.gz variant_mac_maf_sorted.vcf.gz >> bcftools.log 2>&1
                mv variant_mac_maf.recode.vcf.gz FinalSNPs_variant.vcf.gz
                mv invariant.recode.vcf.gz FinalSNPs_invariant.vcf.gz
                rm variant* invariant*
        fi
        
        if [ -e indFiltSNPs.recode.vcf ]; then rm indFiltSNPs.recode.vcf; fi
        rm SNP_temp*
        bgzip -@ $threads $inputVCF
fi

#--- Ending message ---#
echo ""
echo "------------------------------------------------------------------------------------------"
echo "SeqForge vcfpipe is now finished. Thank you for using the program."
echo "For any questions or to report issues, contact Lionel Di Santo at lionel.disanto@unibas.ch"
echo ""
