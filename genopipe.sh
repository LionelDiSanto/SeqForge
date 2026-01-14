#!/bin/bash

#-----------------#
# genopipe v1.0.0 #
#-----------------#

#--- Dependencies ---#
# samtools
# bcftools
# BWA
# fastp
# bowtie2
# seqkit
# htslib
# R

#-----------#
# Code Core #
#-----------#

#--- Program Help and definition of all flags ---#
usage='
genopipe version 1.0.0 Copyright (C) 2025 Lionel Di Santo

genopipe is a bash-based module of SeqForge for the sorting, synchronization, and filtering of raw (demultiplexed) paired-end (PE) read files produced by Illumina technologies with CASAVA 1.8+ FASTQ headers, the mapping of synchronized and filtered reads to a user-provided reference genome, the removal of PCR duplicates from alignment files, and finally the calling of genetic variants (SNPs and Indels).

To run properly, genopipe needs the following dependencies:
- samtools
- bcftools
- BWA
- fastp
- bowtie2
- seqkit
- htslib
- R

** Note that dependencies are not provided with the module and must either be downloaded or loaded as modules by the user **

Before running the genopipe module, make sure that SeqForge and genopipe dependencies are added to your path (if installed locally) using the command lines below:
export PATH=/path_to_dependencies/Dependencies/:$PATH
export PATH=/path_to_SeqForge/SeqForge/:$PATH
!!! If dependencies are loaded as modules, you need not export any paths !!!

genopipe uses specific patterns in file names to run and hence all sequence files **MUST** follow a specific naming convention:

If read files provided are not yet filtered:
PopID_sampID.F.fq.gz (read 1 - forward)
PopID_sampID.R.fq.gz (read 2 - reverse)

If read files provided are already filtered:
PopID_sampID.R1.fq.gz (read 1 - forward)
PopID_sampID.R2.fq.gz (read 2 - reverse)

If files are already aligned (input alignment files must be sorted and in bam format):
PopID_sampID.bam

Note that sequencing files (fastq format) must be compressed using the gzip algorithm.

USAGE:
SeqForge genopipe [options] <reference_genome> <working_directory> <mode>
SeqForge genopipe version (to print the version of the software)

The reference genome (reference_genome) must be placed within the working directory (working_directory) for the module to be able to run read mapping and SNP calling.

As mentioned above, this module can perform read filtering, read mapping, and SNP calling either independently or combined. This feature is determined by the mode.
Six different modes exist:
1/ all - filtering, mapping, and SNPcalling will all be performed sequentially.
2/ filtering - Only filtering of raw reads will be performed.
3/ mapping - Only mapping of filtered reads to the reference genome will be perfomed.
4/ SNPcalling - Only the calling of genetic variants will be performed.
5/ filtering_mapping - Both filtering of raw reads and mapping of filtered reads to the reference genome will be performed.
6/ mapping_SNPcalling - Both mapping of filtered reads to the reference genome and calling of genetic variants will be performed.

Mind that the module will automatically sort and synchronize PE sequences from read 1 and read 2 files (if not already done with a previous run) before starting read filtering
or read mapping. The only instance where sorting and synchronization of reads is not performed is when the module is used solely for SNP calling (when input data is in sorted BAM format).
Additionally, the module will automatically remove PCR duplicates from alignment files by default. Use the appropriate option to disable this behavior if you wish to retain duplicates (see below).

Options to be passed to the module:
    -h, --help                          Show this help page.
    -w, --invariant                     When this argument is specified, variant and invariant sites will be kept in the VCF output 
                                        (if not specified, then only variant sites are returned).
    -c, --discordant                    When this argument is specified, discordant alignments will be kept when using option "--bowtie2".
                                        By default, only concordant alignments are kept during the alignment process.
    -b, --bowtie2                       When this argument is specified, bowtie2 instead of BWA mem will be used for read mapping. While 
                                        BWA mem is a fast and performant aligner, it does not leverage information held by PE sequences. 
                                        Bowtie2, by leveraging this information, may provide more accurate and sensitive alignments, but 
                                        this accuracy and sensitivity come at a price. Computational time when using bowtie2 may be (much) 
                                        higher than when using BWA mem. Note, however, that default settings when using bowtie2 are meant to 
                                        be conservative and may therefore lead to the loss (discard) of many reads during the alignment process. 
                                        When a more lenient mapping is preferred using bowtie2, consider using the -c, --discordant option. 
    -u, --1m-concordant                 By default, when using -b, --bowtie2, all concordant alignments are kept. Specifying this argument enables 
                                        the filtering of all BAM files so that only concordant alignments where reads mapped EXACTLY 1 time to the 
                                        reference genomes are kept (see -q, --min-mapQual for more information). Note that -c, --discordant and -u,
                                        --1m-concordant are MUTALLY EXCLUSIVE.
    -x, --no-dedup                      When this argument is specified, PCR duplicates will NOT be removed from alignment files. By default, PCR 
                                        duplicates are removed.
    -g, --trim-polyG                    If, -g, --trim-polyG is specified, then poly-G tails at the end of reads will be forcefully trimmed. This extra
                                        read filtering step is highly recommended when working with Illumina NextSeq/NovaSeq data.
    -t, --threads           <integer>   Set the number of threads to use for analysis [default: 1].
    -p, --ploidy            <integer>   Ploidy of samples assumed during variant calling [default: 2 (diploid)].
    -d, --depth             <integer>   Number of reads considered per BAM file during variant calling [default: 250].
    -s, --prior             <float>     Expected substitution rate assumed during variant calling. Providing a value of zero (0) disables the use of
                                        this prior during variant calling [default: 0].
    -m, --match             <integer>   Matching score to be passed to aligner [default: 1]. Note that this parameter is needed only when using BWA mem 
                                        as the read mapping program.
    -i, --mismatch          <integer>   Mismatch penalty score to be passed to aligner [default: 4].
    -o, --gapopen           <integer>   Gap open penalty score to be passed to aligner [default: 6].
    -q, --min-mapQual       <integer>   If -u, --1m-concordant is specified, then specifying a value > 0 for this argument will implement an additional
                                        filter where alignments with a mapping quality less than this value will be discarded from BAM files [default: 20].
                                        Specifying a value of 0 disables this filter.                                
'
if [[ $# -eq 0 ]]
    then
        echo ""
        echo "$usage"
        echo ""
        exit 0
fi

PARSED_OPTIONS=$(getopt -o hwcbuxgt:p:d:s:m:i:o:q: -l help,invariant,discordant,bowtie2,1m-concordant,no-dedup,trim-polyG,threads:,ploidy:,depth:,prior:,match:,mismatch:,gapopen:,min-mapQual: -- "$@")
if [[ $? -ne 0 ]]; then echo -e "\nError with parsing arguments\n"; exit 1; fi
eval set -- "$PARSED_OPTIONS"

while true; do
    case "$1" in
        -h|--help) 
            echo "$usage"
            exit
            ;;
        -w|--invariant)
            invariant=activated
            shift
            ;;
        -c|--discordant)
            discordant=activated
            shift
            ;;
        -b|--bowtie2)
            bowtie=activated
            shift
            ;;
        -u|--1m-concordant)
            onlyconc=activated
            shift
            ;;
        -x|--no-dedup)
            nodedup=activated
            shift
            ;;
        -g|--trim-polyG)
            trimG=activated
            shift
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -p|--ploidy)
            ploidy="$2"
            shift 2
            ;;
        -d|--depth)
            depth="$2"
            shift 2
            ;;
        -s|--prior)
            prior="$2"
            shift 2
            ;;
        -m|--match)
            match="$2"
            shift 2
            ;;
        -i|--mismatch)
            mismatch="$2"
            shift 2
            ;;
        -o|--gapopen)
            gapopen="$2"
            shift 2
            ;;
        -q|--min-mapQual)
            minmapQ="$2"
            shift 2
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
        echo "genopipe version 1.0.0 Copyright (C) 2025 Lionel Di Santo"
        echo ""
        exit 0
fi
if [ $# -lt 3 ]
    then
        echo ""
        echo "The reference genome, the working directory, or the mode for the module are not provided."
        echo ""
        exit 1
fi
if [ $# -gt 3 ]
    then
        echo ""
        echo "Too many arguments provided! Three arguments only must be provided following [options]. The reference genome, the working directory, and the mode."
        echo ""
        exit 1
fi
if [ $3 != all ] && [ $3 != filtering ] && [ $3 != mapping ] && [ $3 != SNPcalling ] && [ $3 != filtering_mapping ] && [ $3 != mapping_SNPcalling ]
    then
        echo ""
        echo "The argument for 'mode' is incorrect. See help with SeqForge genopipe -h to see all the different options for 'mode'."
        echo ""
        exit 1
fi

#--- Setting the working directory and mode for the module---#
cd $2
mode=$3

#--- Welcome message ---#
echo ""
echo "#-----------------------------------------------------------#"
echo "# genopipe version 1.0.0 Copyright (C) 2025 Lionel Di Santo #"
echo "#-----------------------------------------------------------#"
program=("samtools" "bcftools" "bwa" "fastp" "htsfile" "bowtie2" "seqkit" "R")
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

#--- Print run parameters ---#
echo ""
echo "Parameters defined for genopipe run"
echo "-----------------------------------"
echo ""
#-- General run information --#
echo "--GENERAL RUN INFORMATION--"
if [ $mode = all ] || [ $mode = mapping ] || [ $mode = filtering_mapping ] || [ $mode = mapping_SNPcalling ]
    then
        if [ -z $bowtie ]; then echo "Aligner used: BWA mem"; fi
        if [ ! -z $bowtie ]; then echo "Aligner used: Bowtie2"; fi
fi
echo "Provided reference genome: $1"
echo "Mode desired for the module: $mode"
echo "Working directory provided: $(pwd)"
if [ -z $threads ]
    then
        echo "Number of threads requested for analysis: 1"
        threads=1
    else
        echo "Number of threads requested for analysis: $threads"
fi
echo ""

#-- Parameters set for raw read filtering --#
if [ $mode = all ] || [ $mode = filtering ] || [ $mode = filtering_mapping ]
    then
    if [ ! -z $trimG ]
        then
            echo "--PARAMETERS SET FOR RAW READ FILTERING--"
            echo "Option -g,--trim-polyG activated: poly-G tails at the end of reads will be removed during read filtering."
    fi
    echo ""
fi

#-- Parameters set for read alignment to the reference genome --#
if [ $mode = all ] || [ $mode = mapping ] || [ $mode = filtering_mapping ] || [ $mode = mapping_SNPcalling ]
    then
    echo "--PARAMETERS SET FOR READ ALIGNMENT TO THE REFERENCE GENOME--"
    if [ ! -z $bowtie ]
        then
            echo "Option -b, --bowtie2 activated: Bowtie2 will be used for read mapping instead of BWA mem."
            aln=Bowtie2
        else
            aln="BWA mem"
    fi
    if [ ! -z $discordant ] && [ ! -z $bowtie ]
        then
            echo "Option -c, --discordant activated: Discordant alignments will be authorized during read mapping with Bowtie2."
    fi
    if [ ! -z $onlyconc ] && [ ! -z $bowtie ]
        then
            if [ -z $minmapQ ]; then minmapQ=20; fi
            echo "Option -u, --1m-concordant activated: Only alignments where reads mapped EXACTLY 1 time will be kept in BAM files following read mapping with Bowtie2."
            if [ $minmapQ -gt 0 ]; then echo "(Additionally, alignments with a mapping quality < $minmapQ will be discarded from all BAM files)"; fi
            if [ ! -z $discordant ]; then echo ""; echo "-c, --discordant CANNOT be used together with -u, --1m-concordant"; echo ""; exit 1; fi
    fi
    if [ ! -z $nodedup ]
        then
            echo "Option -x, --no-dedup activated: Duplicated read pairs and PCR duplicates will NOT be removed from filtered read and alignment files."
    fi
    if [ ! -z $match ] && [ -z $bowtie ]
        then
            echo "Matching score to be passed to aligner: $match"
    fi
    if [ -z $match ] && [ -z $bowtie ]
        then
            echo "Matching score to be passed to aligner: 1"
            match=1
    fi
    if [ -z $mismatch ]
        then
            if [ -z $bowtie ]
                then
                    echo "Mismatch penalty score to be passed to aligner: 4"
                    mismatch=4
            fi
            if [ ! -z $bowtie ]
                then
                    echo "Mismatch penalty (MX) score to be passed to aligner: 4"
                    mismatch=4
                    echo "Mismatch penalty (MN) score to be passed to aligner: $(Rscript -e "cat(round($mismatch/3,0))")"
            fi
        else
            if [ -z $bowtie ]
                then
                    echo "Mismatch penalty score to be passed to aligner: $mismatch"
            fi
            if [ ! -z $bowtie ]
                then
                    echo "Mismatch penalty (MX) score to be passed to aligner: $mismatch"
                    echo "Mismatch penalty (MN) score to be passed to aligner: $(Rscript -e "cat(round($mismatch/3,0))")"
            fi
    fi
    if [ -z $gapopen ]
        then
            if [ -z $bowtie ]
                then  
                    echo "Gap open penalty score to be passed to aligner: 6"
                    gapopen=6
            fi
            if [ ! -z $bowtie ]
                then
                    echo "Gap open penalty score to be passed to aligner: 6"
                    gapopen=6
                    echo "Extension penalty score to be passed to aligner: $(Rscript -e "cat(round($gapopen/1.7,0))")"
            fi
        else
            if [ -z $bowtie ]
                then
                    echo "Gap open penalty score to be passed to aligner: $gapopen"
            fi
            if [ ! -z $bowtie ]
                then
                    echo "Gap open penalty score to be passed to aligner: $gapopen"
                    echo "Extension penalty score to be passed to aligner: $(Rscript -e "cat(round($gapopen/1.7,0))")"
            fi
    fi
    echo ""
fi

#-- Parameters set for genetic variants calling --#
if [ $mode = all ] || [ $mode = SNPcalling ] || [ $mode = mapping_SNPcalling ]
    then
    echo "--PARAMETERS SET FOR GENETIC VARIANTS CALLING--"
    if [ ! -z $invariant ]
        then
            echo "Option -w, --invariant activated: Invariant sites will be kept in the VCF file output."
    fi
    if [ -z $ploidy ]
        then
            echo "Ploidy of samples assumed: diploid (ploidy = 2)"
            ploidy=2
        else
            echo "Ploidy of samples assumed: ploidy = $ploidy"
    fi
    if [ -z $depth ]
        then
            echo "Number of reads considered per BAM file for SNP calling: 250"
            depth=250
        else
            echo "Number of reads considered per BAM file for SNP calling: $depth"
    fi
    if [ -z $prior ]
        then
            echo "No prior expectation on substitution rate assumed during SNP calling (default value is 0)"
            prior=0
        else
            echo "Expected substitution rate assumed during SNP calling: $prior"
    fi
    echo ""
fi

#--- Synchronizing paired-end files using seqkit ---#
echo "Sorting and synchronizing paired-end files using seqkit"
echo "-------------------------------------------------------"
if [ $(ls . | grep -E "paired.*\.fq\.gz" | wc -l) -eq 0 ] && [ ! $mode = SNPcalling ]
    then
        if [ $(ls . | grep -E ".*F\.fq\.gz" | wc -l) -gt 0 ]; then nf=$(ls *F.fq.gz | wc -l); else nf=0; fi
        if [ $(ls . | grep -E ".*R\.fq\.gz" | wc -l) -gt 0 ]; then nr=$(ls *R.fq.gz | wc -l); else nr=0; fi
        if [ $(ls . | grep -E ".*R1\.fq\.gz" | wc -l) -gt 0 ]; then n1=$(ls *R1.fq.gz | wc -l); else n1=0; fi
        if [ $(ls . | grep -E ".*R2\.fq\.gz" | wc -l) -gt 0 ]; then n2=$(ls *R2.fq.gz | wc -l); else n2=0; fi
        if [ $nf -eq 0 ] && [ $nr -eq 0 ] && [ $n1 -eq 0 ] && [ $n2 -eq 0 ]
            then
                echo ""
                echo "File with names F/R1 and R/R2 (representing raw or filtered read files) are absent. Please check the data provided and the naming convention (see SeqForge genopipe -h for more details)."
                echo ""
                exit 1
        fi
        if [ $nf -gt 0 ] && [ $nr -gt 0 ] && [ $n1 -eq 0 ] && [ $n2 -eq 0 ]
            then
                b1="F"
                b2="R"
        fi
        if [ $nf -eq 0 ] && [ $nr -eq 0 ] && [ $n1 -gt 0 ] && [ $n2 -gt 0 ]
            then
                b1="R1"
                b2="R2"
        fi
        if [ $nf -gt 0 ] && [ $nr -gt 0 ] && [ $n1 -gt 0 ] && [ $n2 -gt 0 ]
            then
                b1="R1"
                b2="R2"
        fi
        if [ ! $(ls *${b1}.fq.gz | wc -l) -eq $(ls *${b2}.fq.gz | wc -l) ]
            then
                echo ""
                echo "All files do not appear to have their corresponding pair. Please check the data provided."
                echo ""
                exit 1
        fi
        echo "Number of individuals detected:" $(ls *${b1}.fq.gz | wc -l) "(2 times this number of sequence files)."
        echo ""
        echo "Checking for sorted and synchronized files"
        ls *${b1}.fq.gz | sed "s/.${b1}.fq.gz//g" > ind.list
        for line in $(cat ind.list)
            do
                zcat $line.${b1}.fq.gz | grep "@.*[YN]" | sed -E 's/([12])(:[NY])/x\2/g' > Fnames & 
                zcat $line.${b2}.fq.gz | grep "@.*[YN]" | sed -E 's/([12])(:[NY])/x\2/g' > Rnames &
                wait
                diff Fnames Rnames > diff.names
                if [ $(cat diff.names | wc -l) -eq 0 ]
                    then 
                        echo -e "\tReads for individual $line are synchronized. Adding prefix 'paired' in front of read files."
                        mv $line.${b1}.fq.gz paired_$line.${b1}.fq.gz
                        mv $line.${b2}.fq.gz paired_$line.${b2}.fq.gz
                    else
                        echo "\tReads for individual $line are NOT synchronized. Sorting and synchronization will be performed for this individual."
                        echo -e "$line.${b1}.fq.gz\n$line.${b2}.fq.gz" >> unpaired.list
                fi
                rm Fnames Rnames diff.names
        done
        rm ind.list
        if [ -e unpaired.list ]
            then
                echo ""
                echo "Sorting and synchronizing paired-end files using seqkit will performed for the following files:"
                echo "$(cat unpaired.list | tr '\n' ',' | sed 's/,$//')."
                echo ""
                echo "Sorting FORWARD (F/R1) read files"
                cat unpaired.list | grep "${b1}.fq.gz" > raw.list
                while [ $(cat raw.list | wc -l) -gt 0 ]
                    do
                        cat raw.list | head -n $threads > temp.list
                        count=0
                        while read -r line
                            do
                                count=$(($count+1))
                                if [ $count -eq $(cat temp.list | wc -l) ]
                                    then
                                        seqkit sort -n $line -o sorted_$line -j 1 --quiet
                                        wait
                                    else
                                        seqkit sort -n $line -o sorted_$line -j 1 --quiet &
                                fi
                        done < temp.list
                        rm temp.list
                        cat raw.list | tail -n +$(($threads+1)) > raw.list.tmp && mv raw.list.tmp raw.list  
                done
                rm raw.list
                echo "Sorting REVERSE (R/R2) read files"
                cat unpaired.list | grep "${b2}.fq.gz" > raw.list
                while [ $(cat raw.list | wc -l) -gt 0 ]
                    do
                        cat raw.list | head -n $threads > temp.list
                        count=0
                        while read -r line
                            do
                                count=$(($count+1))
                                if [ $count -eq $(cat temp.list | wc -l) ]
                                    then
                                        seqkit sort -n $line -o sorted_$line -j 1 --quiet
                                        wait
                                    else
                                        seqkit sort -n $line -o sorted_$line -j 1 --quiet &
                                fi
                        done < temp.list
                        rm temp.list
                        cat raw.list | tail -n +$(($threads+1)) > raw.list.tmp && mv raw.list.tmp raw.list  
                done
                rm raw.list unpaired.list
                echo "Synchronization of FORWARD and REVERSE read files"
                ls sorted*${b1}.fq.gz > raw.list
                while [ $(cat raw.list | wc -l) -gt 0 ]
                    do
                        cat raw.list | head -n $threads > temp.list
                        count=0
                        while read -r line
                            do
                                count=$(($count+1))
                                if [ $count -eq $(cat temp.list | wc -l) ]
                                    then
                                        seqkit pair -1 $line -2 $(basename $line | sed "s/${b1}.fq.gz/${b2}.fq.gz/") -j 1 --quiet
                                        wait
                                    else
                                        seqkit pair -1 $line -2 $(basename $line | sed "s/${b1}.fq.gz/${b2}.fq.gz/") -j 1 --quiet &
                                fi
                        done < temp.list
                        rm temp.list
                        cat raw.list | tail -n +$(($threads+1)) > raw.list.tmp && mv raw.list.tmp raw.list  
                done
                rm raw.list
                ls . | grep -v -E "^sorted|input|$(basename $1)" > input
                mkdir input_reads
                while read -r line; do mv $line input_reads/$line; done < input
                rm input sorted*.${b1}.fq.gz sorted*.${b2}.fq.gz
                echo ""
                echo "Renaming files to add prefix 'paired'"
                for file in sorted*.fq.gz; do mv $file $(basename $file | sed 's/.paired//'); done
                for file in sorted*.fq.gz; do  echo "$(basename $file | sed 's/sorted_//') --> $(basename $file | sed 's/sorted_/paired_/')"; mv $file $(basename $file | sed 's/sorted_/paired_/'); done
        fi
    else
        if [ ! $mode = SNPcalling ]; then echo "Synchronized files detected ($(ls paired*${b1}.fq.gz | wc -l)) - SKIPPING sorting and sychronization of paired-end read files."; fi
        if [ $mode = SNPcalling ]; then echo "Sorting and synchronization of read files is not needed anymore when running the module only in SNPcalling mode - SKIPPING sorting and sychronization of paired-end read files."; fi
fi

#--- Raw read filtering using fastp ---#
if [ $mode = all ] || [ $mode = filtering ] || [ $mode = filtering_mapping ]
then
    echo ""
    echo "Raw reads filtering using fastp"
    echo "-------------------------------"
    ls *F.fq.gz > raw.list
    while [ $(cat raw.list | wc -l) -gt 0 ]
        do
            cat raw.list | head -n $threads > temp.list
            echo -e "Processing individuals:\n $(cat temp.list | sed 's/.F.fq.gz//g' | sed -e "s/^/\t/g")"
            count=0
            while read -r line
                do
                    count=$(($count+1))
                    if [ $count -eq $(cat temp.list | wc -l) ]
                        then
                            if [ -z $trimG ]
                                then 
                                    fastp --in1 $line --in2 $(basename $line | sed 's/F.fq.gz/R.fq.gz/') --out1 $(basename $line | sed 's/F.fq.gz/R1.fq.gz/') --out2 $(basename $line | sed 's/F.fq.gz/R2.fq.gz/') --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 --correction -q 15 -u 50 -h $(basename $line | sed 's/F.fq.gz/html/') -j $(basename $line | sed 's/F.fq.gz/json/') --detect_adapter_for_pe > $(basename $line | sed 's/F.fq.gz/log/') 2>&1
                                    wait
                                else
                                    fastp --in1 $line --in2 $(basename $line | sed 's/F.fq.gz/R.fq.gz/') --out1 $(basename $line | sed 's/F.fq.gz/R1.fq.gz/') --out2 $(basename $line | sed 's/F.fq.gz/R2.fq.gz/') --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 --correction --trim_poly_g -q 15 -u 50 -h $(basename $line | sed 's/F.fq.gz/html/') -j $(basename $line | sed 's/F.fq.gz/json/') --detect_adapter_for_pe > $(basename $line | sed 's/F.fq.gz/log/') 2>&1
                                    wait
                            fi
                        else
                            if [ -z $trimG ]
                                then
                                    fastp --in1 $line --in2 $(basename $line | sed 's/F.fq.gz/R.fq.gz/') --out1 $(basename $line | sed 's/F.fq.gz/R1.fq.gz/') --out2 $(basename $line | sed 's/F.fq.gz/R2.fq.gz/') --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 --correction -q 15 -u 50 -h $(basename $line | sed 's/F.fq.gz/html/') -j $(basename $line | sed 's/F.fq.gz/json/') --detect_adapter_for_pe > $(basename $line | sed 's/F.fq.gz/log/') 2>&1 &
                                else
                                    fastp --in1 $line --in2 $(basename $line | sed 's/F.fq.gz/R.fq.gz/') --out1 $(basename $line | sed 's/F.fq.gz/R1.fq.gz/') --out2 $(basename $line | sed 's/F.fq.gz/R2.fq.gz/') --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 --correction --trim_poly_g -q 15 -u 50 -h $(basename $line | sed 's/F.fq.gz/html/') -j $(basename $line | sed 's/F.fq.gz/json/') --detect_adapter_for_pe > $(basename $line | sed 's/F.fq.gz/log/') 2>&1 &
                            fi
                    fi
            done < temp.list
            rm temp.list
            cat raw.list | tail -n +$(($threads+1)) > raw.list.tmp && mv raw.list.tmp raw.list
            if [ $(cat raw.list | wc -l) -gt 0 ]; then echo ""; fi
    done
    rm raw.list
    mkdir fastp_trim_reports fastp_logfiles
    mv *.html *.json fastp_trim_reports
    mv *.log fastp_logfiles
fi

#--- Aligning reads to the reference genome ---#
if [ $mode = all ] || [ $mode = mapping ] || [ $mode = filtering_mapping ] || [ $mode = mapping_SNPcalling ]
then
    echo ""
    echo "Reads alignment to the reference genome using $aln"
    echo "-----------------------------------------------------"
    n1=$(ls . | grep -E ".*R1\.fq\.gz" | wc -l)
    n2=$(ls . | grep -E ".*R2\.fq\.gz" | wc -l)
    if [ $n1 -eq 0 ] || [ $n2 -eq 0 ]
        then
            echo ""
            echo "File with names R1 and R2 (representing filtered read files) are absent. Either run the module in filtering mode to generate these files or rename your read files according to the naming convention (see SeqForge genopipe -h for more details)."
            echo ""
            exit 1
    fi
    if [ ! $n1 -eq $n2 ]
        then
            echo ""
            echo "All files do not appear to have their corresponding pair. Please check the data provided."
            echo ""
            exit 1
    fi
    if [ ! -e $1 ]
        then
            echo ""
            echo "The reference genome is absent from the working directory."
            echo ""
            exit 1
    fi
    if [ -z $bowtie ]
        then
            if [ $(ls . | grep "$1.bwt" | wc -l) -eq 0 ]; then echo "Indexing of the reference genome: $1"; bwa index $1; echo ""; fi
            if [ $(ls . | grep "$1.fai" | wc -l) -eq 0 ]; then samtools faidx $1; fi
            c=0
            for file in *R1.fq.gz
                do
                    if [ $c -eq 0 ]; then echo "Processing individual: $(basename $file | sed 's/.R1.fq.gz//')"; else echo -e "\nProcessing individual: $(basename $file | sed 's/.R1.fq.gz//')"; fi
                    rgid=$(zcat $file | head -n 1 | awk -v sname=$(basename $file | sed 's/.R1.fq.gz//' | sed 's/paired_//') -F':' '{print "@RG\\tID:"$3"."$4"\\tSM:"sname"\\tPL:ILLUMINA\\tPU:"$3"."$4"."$NF}')
                    bwa mem -t $threads -v 2 -A $match -B $mismatch -O $gapopen -R $rgid $1 $file $(basename $file | sed 's/R1.fq.gz/R2.fq.gz/') > $(basename $file | sed 's/R1.fq.gz/sam/')
                    sed -i '/^$/d' $(basename $file | sed 's/R1.fq.gz/sam/')
                    samtools view -h -b --fast -t $1.fai -S -o $(basename $file | sed 's/R1.fq.gz/bam.bai/') $(basename $file | sed 's/R1.fq.gz/sam/') 
                    samtools sort -@ $threads -o $(basename $file | sed 's/R1.fq.gz/bam/') $(basename $file | sed 's/R1.fq.gz/bam.bai/') 
                    ((c++))
            done
            rm *.sam
            rm *.bam.bai
    fi
    if [ ! -z $bowtie ]
        then
            if [ $(ls . | grep -E "$1.*\.bt2" | wc -l) -eq 0 ]; then echo -e "Indexing of the reference genome: $1 \n"; bowtie2-build -f --quiet --threads $threads $1 $1; fi
            if [ $(ls . | grep "$1.fai" | wc -l) -eq 0 ]; then samtools faidx $1; fi
            c=0
            for file in *R1.fq.gz
                do
                    if [ $c -eq 0 ]; then echo "Processing individual: $(basename $file | sed 's/.R1.fq.gz//')"; echo " $(basename $file | sed 's/.R1.fq.gz//')" >> bowtie2.log; else echo -e "\nProcessing individual: $(basename $file | sed 's/.R1.fq.gz//')"; echo -e "\n $(basename $file | sed 's/.R1.fq.gz//')" >> bowtie2.log; fi
                    if [ -z $discordant ]; then bowtie2 -p $threads -x $1 -1 $file -2 $(basename $file | sed 's/R1.fq.gz/R2.fq.gz/') -S $(basename $file | sed 's/R1.fq.gz/sam/') --rg-id $(zcat $file | head -n 1 | awk -F':' '{print $3"."$4}') --rg "SM:$(basename $file | sed 's/.R1.fq.gz//' | sed 's/paired_//')" --rg "PL:ILLUMINA" --rg $(zcat $file | head -n 1 | awk -F':' '{print "PU:"$3"."$4"."$NF}') -q --very-sensitive --no-1mm-upfront --end-to-end --mp $mismatch,$(Rscript -e "cat(round($mismatch/3,0))") --np 1 --rdg $gapopen,$(Rscript -e "cat(round($gapopen/1.7,0))") --rfg $gapopen,$(Rscript -e "cat(round($gapopen/1.7,0))") --fr --no-mixed --no-contain --no-discordant --no-unal >> bowtie2.log 2>&1 ; fi
                    if [ ! -z $discordant ]; then bowtie2 -p $threads -x $1 -1 $file -2 $(basename $file | sed 's/R1.fq.gz/R2.fq.gz/') -S $(basename $file | sed 's/R1.fq.gz/sam/') --rg-id $(zcat $file | head -n 1 | awk -F':' '{print $3"."$4}') --rg "SM:$(basename $file | sed 's/.R1.fq.gz//' | sed 's/paired_//')" --rg "PL:ILLUMINA" --rg $(zcat $file | head -n 1 | awk -F':' '{print "PU:"$3"."$4"."$NF}') -q --very-sensitive --no-1mm-upfront --end-to-end --mp $mismatch,$(Rscript -e "cat(round($mismatch/3,0))") --np 1 --rdg $gapopen,$(Rscript -e "cat(round($gapopen/1.7,0))") --rfg $gapopen,$(Rscript -e "cat(round($gapopen/1.7,0))") --fr --no-mixed --no-contain --no-unal >> bowtie2.log 2>&1; fi
                    samtools view -h -b --fast -t $1.fai -S -o $(basename $file | sed 's/R1.fq.gz/bam.bai/') $(basename $file | sed 's/R1.fq.gz/sam/') 
                    samtools sort -@ $threads -o $(basename $file | sed 's/R1.fq.gz/bam/') $(basename $file | sed 's/R1.fq.gz/bam.bai/') 
                    echo "" >> bowtie2.log
                    ((c++))
            done
            rm *.sam
            rm *.bam.bai
            if [ ! -z $onlyconc ] && [ -z $discordant ]
                then
                    echo ""
                    if [ $minmapQ -eq 0 ]
                        then 
                            echo "Filtering of BAM files to keep uniquely mapped concordant alignments"
                            echo "--------------------------------------------------------------------"
                        else
                            echo "Filtering of BAM files to keep uniquely mapped concordant alignments with a mapping quality >= $minmapQ"
                            echo "-------------------------------------------------------------------------------------------------------"
                    fi
                    c=0
                    for file in *bam
                        do
                            if [ $c -eq 0 ]; then echo "Processing: $file"; else echo -e "\nProcessing: $file"; fi
                            samtools view -h -F 256 -F 2048 -q $minmapQ -o $(basename $file | sed 's/bam/bam.conc/') $file 
                            echo "  Kept $(samtools view -c $(basename $file | sed 's/bam/bam.conc/')) of $(samtools view -c $file) alignments."
                            ((c++))
                    done
                    rm *bam
                    for file in *bam.conc; do mv $file $(basename $file | sed 's/.conc//'); done
            fi
    fi
    echo ""
    echo "Removing PCR duplicates from alignment files (BAM files) using samtools"
    echo "-----------------------------------------------------------------------"
    if [ -z $nodedup ]
       then
            nb=$(ls . | grep ".bam" | wc -l)
            if [ $nb -eq 0 ]
                then
                    echo ""
                    echo "No bam files are detected within the working directory. If these files are absent, think of running the module in mapping mode."
                    echo ""
                    exit 1
            fi
            ls *.bam > mybamlist.txt
            while read -r line
                do
                    echo "Processing $line --> $(basename $line | sed 's/bam/markdup.bam/')"
                    samtools collate -@ $threads -o $(basename $line | sed 's/bam/namecollate.bam/') $line 
                    samtools fixmate -@ $threads -m $(basename $line | sed 's/bam/namecollate.bam/') $(basename $line | sed 's/bam/fixmate.bam/') 
                    samtools sort -@ $threads -o $(basename $line | sed 's/bam/positionsort.bam/') $(basename $line | sed 's/bam/fixmate.bam/') 
                    samtools markdup -@ $threads -d 0 -r -c -f $(basename $line | sed 's/bam/dupstats/') $(basename $line | sed 's/bam/positionsort.bam/') $(basename $line | sed 's/bam/markdup.bam/')
                    echo "" 
            done < mybamlist.txt
            rm *namecollate.bam *fixmate.bam *positionsort.bam mybamlist.txt
            mkdir Duplicated DeDupStats
            mv *bam Duplicated
            mv *dupstats DeDupStats
            mv ./Duplicated/*markdup.bam .
        else
            echo "-x, --no-dedup activated: Skipping deduplication of alignment files with samtools."
            echo ""
    fi
fi

#--- SNPs calling ---#
if [ $mode = all ] || [ $mode = SNPcalling ] || [ $mode = mapping_SNPcalling ]
then
    if [ $mode = SNPcalling ]; then echo ""; fi
    echo "SNP calling using bcftools"
    echo "--------------------------"
    nb=$(ls . | grep ".bam" | wc -l)
    if [ $nb -eq 0 ]
        then
            echo ""
            echo "No bam files are detected within the working directory. If these files are absent, think of running the module in mapping mode."
            echo ""
            exit 1
    fi
    if [ ! -e $1 ]
       then
           echo ""
           echo "The reference genome is absent from the working directory."
           echo ""
           exit 1
    fi
    ls *.bam > mybamlist.txt
    for i in ./*.bam; do samtools view -H $i | grep '@RG' | cut -f3 | cut -f2 -d: >> join1; done
    cat join1 | sed 's/_[0-9]*$//' >> join2
    paste join1 join2 > groups.txt
    rm join1 join2
    if [ -z $invariant ]
        then
            bcftools mpileup -f $1 --threads $threads --max-depth $depth --bam-list ./mybamlist.txt -a DP,QS,AD --full-BAQ --config illumina | bcftools call -o RawSNPs.vcf.gz -Oz --threads $threads --ploidy $ploidy --variants-only -m --prior $prior -G groups.txt
        else
            bcftools mpileup -f $1 --threads $threads --max-depth $depth --bam-list ./mybamlist.txt -a DP,QS,AD --full-BAQ --config illumina | bcftools call -o RawSNPs.vcf.gz -Oz --threads $threads --ploidy $ploidy -m --prior $prior -G groups.txt
    fi
fi

#--- Ending message ---#
echo ""
echo "------------------------------------------------------------------------------------------"
echo "SeqForge genopipe is now finished. Thank you for using the program."
echo "For any questions or to report issues, contact Lionel Di Santo at lionel.disanto@unibas.ch"
echo ""
