#!/bin/bash

echo "#---------------------------------------------------------------------"
echo "# This script is adopted from run-gen-ref in bwa.kit 0.7.12 by Heng Li"
echo "# available from https://github.com/lh3/bwa/tree/master/bwakit"
echo "#"
echo "# run this to download GRCh38 reference."
echo "#"
echo "#---------------------------------------------------------------------"
echo

pushd `dirname $0` > /dev/null
SCRIPTD=`pwd`
popd > /dev/null
resources_dir=$SCRIPTD/../resources

url38NoAltDecoy="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
url38d="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz"
url38="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"



function print_usage {
    echo "Usage: $0 <hs38|hs38D|hs38DH|hs38NoAltDH> "
    echo ""
    echo "      The file containting the HLA sequences is located at resources/hla_bwa.kit.fna.gz"
    echo ""
    echo "Analysis sets:"
    echo "  hs38         primary assembly of GRCh38 (incl. chromosomes, unplaced and unlocalized contigs) and EBV"
    echo "  hs38D        hs38 + decoy contigs"
    echo "  hs38DH       hs38 + ALT contigs + decoy contigs + HLA genes (recommended for GRCh38 mapping)"
    echo "  hs38NoAltDH  hs38 + decoy contigs + HLA alleles from [bwa.kit] (recommended for HLA typing)"
    echo ""
    echo "Note: This script downloads human reference genomes from NCBI ftp server"
    exit 1;
}
#if [ $# -eq 2 ]; then
#    if [ ! -e $2 ]; then
#	echo "ERROR: $2 NOT FOUND!"
#	print_usage
#    fi
#fi

if [ $# -eq 1 ]; then
    if [ $1 == "hs38DH" ]; then
	if [ ! -e $resources_dir/$1.fa ]; then
	    (wget -O- $url38d | gzip -dc; gzip -dc $resources_dir/hla_bwa.kit.fna.gz) > $resources_dir/$1.fa
	else
	    print_usage
	fi
    elif [ $1 == "hs38" ]; then
	if [ $# -eq 1 ]; then
	    (wget -O- $url38 | gzip -dc) > $resources_dir/$1.fa
	else
	    print_usage
	fi
    elif [ $1 == "hs38D" ]; then
	if [ $# -eq 1 ]; then
	    (wget -O- $url38d | gzip -dc) > $resources_dir/$1.fa
	else
	    print_usage
	fi
    elif [ $1 == "hs38NoAltDH" ]; then
	if [ $# -eq 1 ]; then
	    (wget -O- $url38NoAltDecoy | gzip -dc; gzip -dc $resources_dir/hla_bwa.kit.fna.gz) > $resources_dir/$1.fa
	else
	    print_usage
	fi
    else
	echo "ERROR: unknown genome build"
	echo
	print_usage
    fi
    
    [ ! -f $resources_dir/$1.fa.bwt ] && echo -e "\nPlease run 'bwa index $resources_dir/$1.fa'...\n"
else
    print_usage
fi
