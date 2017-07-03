#
# Part of Kourami HLA typer/assembler
# (c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
# See LICENSE for licensing.
#

#!/bin/bash

pushd `dirname $0` > /dev/null
SCRIPTD=`pwd`
popd > /dev/null

samtools_sort_memory_per_thread=2G
num_processors=8
kourami_db=$SCRIPTD/../db
me=`basename $0`


function usage {
    echo "HLA-related reads extractor for Kourami"
    echo "Note: Use this if you have bam file aligned to GRCh38 [NoAlt] (primary assembly + decoy + HLA from [bwa-kit] )"
    echo "USAGE: <PATH-TO>/$me -d [Kourami panel db] -r [refGenome] <sample_id> <bamfile>"
    echo
    echo " sample_id        : desired sample name (ex: NA12878) [required]"
    echo
    echo " bamfile          : sorted and indexed bam to hs38NoAltDH (ex: NA12878.bam) [required]"
    echo
    echo "------------------ Optional Parameters -----------------"
    echo " -d [panel DB]    : Path to Kourami panel db. [Default: db directory under Kourami installation kourami/db]"
    echo
    echo " -r [Ref Gemome]  : path to hs38NoAltDH (primary assembly + decoy + HLA [bwa-kit])" 
    echo "                    USE download_grch38.sh script to obtain the reference."
    echo "                    MUST BE BWA INDEXED prior to running this script."
    echo "                    If not given, it assumes, hs38NoAltDH.fa is in resources dir."
    echo
    echo " -h               : print this message."
    echo
    exit 1
}

# print usage when no argument is given
if [ $# -lt 1 ]; then
    usage
fi

while getopts :d:r:h FLAG; do
    case $FLAG in 
	d) 
	    kourami_db=$OPTARG
	    ;;
	h)
	    usage
	    ;;
	\?)
	    echo "Unrecognized option -$OPTARG. See usage:"
	    usage
	    ;;
    esac
done

shift $((OPTIND-1))

if [ $# -lt 2 ]; then
    echo "Missing one or more required arguments."
    usage
fi

sampleid=$1
bam_path=$2

merged_hla_panel=$kourami_db/All_FINAL_with_Decoy.fa.gz
bam_for_kourami=$sampleid\_on_KouramiPanel.bam
samtools_bin=`(which samtools)`
bwa_bin=`(which bwa)`
bamUtil=`(which bam)`
if [ -z "$bamUtil" ]; then
    echo "missing bamUtil";
    echo "bamUtil available from https://github.com/statgen/bamUtil"
    exit 1;
fi
#bamUtil=$HOME/bamUtil_1.0.13/bamUtil-master/bin/bam 

if [ ! -x "$samtools_bin" ] || [ ! -x "$bwa_bin" ] || [ ! -x "$bamUtil" ];then
    echo "Please make sure samtools, bwa, and bamUtil are installed"
    exit 1
fi

if [ ! -e "$bam_path" ] || [ ! -e "$kourami_db" ] || [ ! -e "$merged_hla_panel" ];then
    echo "Missing one of the following files/directories (38DH):\n"
    echo "$bam_path"
    echo "$kourami_db"
    echo "$merged_hla_panel"
    exit 1
fi

echo ">>>>>>>>>>>>>>>> extracting reads mapping to HLA loci and ALT contigs (38DH_NoAlt)"
$samtools_bin view -b $bam_path chr6:29723340-29727296 chr6:29726601-29749049 chr6:29826979-29831122 \
    chr6:29887760-29891080 chr6:29942470-29945884 chr6:30005971-30009956 chr6:30259562-30266951 \
    chr6:30489406-30494205 chr6:31268749-31272136 chr6:31353866-31357245 chr6:31399784-31415316 \
    chr6:31494881-31511124 chr6:32439842-32445051 chr6:32517377-32530229 chr6:32552713-32560002 \
    chr6:32578770-32589836 chr6:32637406-32643652 chr6:32659464-32666689 chr6:32741386-32746887 \
    chr6:32756098-32763553 chr6:32812763-32817048 chr6:32821833-32838770 chr6:32845209-32852787 \
    chr6:32934629-32941070 chr6:32948614-32953122 chr6:33004183-33009612 chr6:33064569-33080778 \
    chr6:33075926-33089696 chr6:33112516-33129113 HLA-A*01:01:01:01: HLA-A*01:01:01:02N: HLA-A*01:01:38L: HLA-A*01:02: HLA-A*01:03: HLA-A*01:04N: HLA-A*01:09: HLA-A*01:11N: HLA-A*01:14: HLA-A*01:16N: HLA-A*01:20: HLA-A*02:01:01:01: HLA-A*02:01:01:02L: HLA-A*02:01:01:03: HLA-A*02:01:01:04: HLA-A*02:02:01: HLA-A*02:03:01: HLA-A*02:03:03: HLA-A*02:05:01: HLA-A*02:06:01: HLA-A*02:07:01: HLA-A*02:10: HLA-A*02:251: HLA-A*02:259: HLA-A*02:264: HLA-A*02:265: HLA-A*02:266: HLA-A*02:269: HLA-A*02:279: HLA-A*02:32N: HLA-A*02:376: HLA-A*02:43N: HLA-A*02:455: HLA-A*02:48: HLA-A*02:51: HLA-A*02:533: HLA-A*02:53N: HLA-A*02:57: HLA-A*02:60:01: HLA-A*02:65: HLA-A*02:68: HLA-A*02:77: HLA-A*02:81: HLA-A*02:89: HLA-A*02:95: HLA-A*03:01:01:01: HLA-A*03:01:01:02N: HLA-A*03:01:01:03: HLA-A*03:02:01: HLA-A*03:11N: HLA-A*03:21N: HLA-A*03:36N: HLA-A*11:01:01: HLA-A*11:01:18: HLA-A*11:02:01: HLA-A*11:05: HLA-A*11:110: HLA-A*11:25: HLA-A*11:50Q: HLA-A*11:60: HLA-A*11:69N: HLA-A*11:74: HLA-A*11:75: HLA-A*11:77: HLA-A*23:01:01: HLA-A*23:09: HLA-A*23:38N: HLA-A*24:02:01:01: HLA-A*24:02:01:02L: HLA-A*24:02:01:03: HLA-A*24:02:03Q: HLA-A*24:02:10: HLA-A*24:03:01: HLA-A*24:07:01: HLA-A*24:08: HLA-A*24:09N: HLA-A*24:10:01: HLA-A*24:11N: HLA-A*24:152: HLA-A*24:20: HLA-A*24:215: HLA-A*24:61: HLA-A*24:86N: HLA-A*25:01:01: HLA-A*26:01:01: HLA-A*26:11N: HLA-A*26:15: HLA-A*26:50: HLA-A*29:01:01:01: HLA-A*29:01:01:02N: HLA-A*29:02:01:01: HLA-A*29:02:01:02: HLA-A*29:46: HLA-A*30:01:01: HLA-A*30:02:01:01: HLA-A*30:02:01:02: HLA-A*30:04:01: HLA-A*30:89: HLA-A*31:01:02: HLA-A*31:01:23: HLA-A*31:04: HLA-A*31:14N: HLA-A*31:46: HLA-A*32:01:01: HLA-A*32:06: HLA-A*33:01:01: HLA-A*33:03:01: HLA-A*33:07: HLA-A*34:01:01: HLA-A*34:02:01: HLA-A*36:01: HLA-A*43:01: HLA-A*66:01:01: HLA-A*66:17: HLA-A*68:01:01:01: HLA-A*68:01:01:02: HLA-A*68:01:02:01: HLA-A*68:01:02:02: HLA-A*68:02:01:01: HLA-A*68:02:01:02: HLA-A*68:02:01:03: HLA-A*68:02:02: HLA-A*68:03:01: HLA-A*68:08:01: HLA-A*68:113: HLA-A*68:17: HLA-A*68:18N: HLA-A*68:22: HLA-A*68:71: HLA-A*69:01: HLA-A*74:01: HLA-A*74:02:01:01: HLA-A*74:02:01:02: HLA-A*80:01:01:01: HLA-A*80:01:01:02: HLA-B*07:02:01: HLA-B*07:05:01: HLA-B*07:06: HLA-B*07:156: HLA-B*07:33:01: HLA-B*07:41: HLA-B*07:44: HLA-B*07:50: HLA-B*08:01:01: HLA-B*08:08N: HLA-B*08:132: HLA-B*08:134: HLA-B*08:19N: HLA-B*08:20: HLA-B*08:33: HLA-B*08:79: HLA-B*13:01:01: HLA-B*13:02:01: HLA-B*13:02:03: HLA-B*13:02:09: HLA-B*13:08: HLA-B*13:15: HLA-B*13:25: HLA-B*14:01:01: HLA-B*14:02:01: HLA-B*14:07N: HLA-B*15:01:01:01: HLA-B*15:01:01:02N: HLA-B*15:01:01:03: HLA-B*15:02:01: HLA-B*15:03:01: HLA-B*15:04:01: HLA-B*15:07:01: HLA-B*15:108: HLA-B*15:10:01: HLA-B*15:11:01: HLA-B*15:13:01: HLA-B*15:16:01: HLA-B*15:17:01:01: HLA-B*15:17:01:02: HLA-B*15:18:01: HLA-B*15:220: HLA-B*15:25:01: HLA-B*15:27:01: HLA-B*15:32:01: HLA-B*15:42: HLA-B*15:58: HLA-B*15:66: HLA-B*15:77: HLA-B*15:83: HLA-B*18:01:01:01: HLA-B*18:01:01:02: HLA-B*18:02: HLA-B*18:03: HLA-B*18:17N: HLA-B*18:26: HLA-B*18:94N: HLA-B*27:04:01: HLA-B*27:05:02: HLA-B*27:05:18: HLA-B*27:06: HLA-B*27:07:01: HLA-B*27:131: HLA-B*27:24: HLA-B*27:25: HLA-B*27:32: HLA-B*35:01:01:01: HLA-B*35:01:01:02: HLA-B*35:01:22: HLA-B*35:02:01: HLA-B*35:03:01: HLA-B*35:05:01: HLA-B*35:08:01: HLA-B*35:14:02: HLA-B*35:241: HLA-B*35:41: HLA-B*37:01:01: HLA-B*37:01:05: HLA-B*38:01:01: HLA-B*38:02:01: HLA-B*38:14: HLA-B*39:01:01:01: HLA-B*39:01:01:02L: HLA-B*39:01:01:03: HLA-B*39:01:03: HLA-B*39:01:16: HLA-B*39:01:21: HLA-B*39:05:01: HLA-B*39:06:02: HLA-B*39:10:01: HLA-B*39:13:02: HLA-B*39:14: HLA-B*39:34: HLA-B*39:38Q: HLA-B*40:01:01: HLA-B*40:01:02: HLA-B*40:02:01: HLA-B*40:03: HLA-B*40:06:01:01: HLA-B*40:06:01:02: HLA-B*40:10:01: HLA-B*40:150: HLA-B*40:40: HLA-B*40:72:01: HLA-B*40:79: HLA-B*41:01:01: HLA-B*41:02:01: HLA-B*42:01:01: HLA-B*42:02: HLA-B*42:08: HLA-B*44:02:01:01: HLA-B*44:02:01:02S: HLA-B*44:02:01:03: HLA-B*44:02:17: HLA-B*44:02:27: HLA-B*44:03:01: HLA-B*44:03:02: HLA-B*44:04: HLA-B*44:09: HLA-B*44:138Q: HLA-B*44:150: HLA-B*44:23N: HLA-B*44:26: HLA-B*44:46: HLA-B*44:49: HLA-B*44:56N: HLA-B*45:01:01: HLA-B*45:04: HLA-B*46:01:01: HLA-B*46:01:05: HLA-B*47:01:01:01: HLA-B*47:01:01:02: HLA-B*48:01:01: HLA-B*48:03:01: HLA-B*48:04: HLA-B*48:08: HLA-B*49:01:01: HLA-B*49:32: HLA-B*50:01:01: HLA-B*51:01:01: HLA-B*51:01:02: HLA-B*51:02:01: HLA-B*51:07:01: HLA-B*51:42: HLA-B*52:01:01:01: HLA-B*52:01:01:02: HLA-B*52:01:01:03: HLA-B*52:01:02: HLA-B*53:01:01: HLA-B*53:11: HLA-B*54:01:01: HLA-B*54:18: HLA-B*55:01:01: HLA-B*55:01:03: HLA-B*55:02:01: HLA-B*55:12: HLA-B*55:24: HLA-B*55:48: HLA-B*56:01:01: HLA-B*56:03: HLA-B*56:04: HLA-B*57:01:01: HLA-B*57:03:01: HLA-B*57:06: HLA-B*57:11: HLA-B*57:29: HLA-B*58:01:01: HLA-B*58:31N: HLA-B*59:01:01:01: HLA-B*59:01:01:02: HLA-B*67:01:01: HLA-B*67:01:02: HLA-B*67:02: HLA-B*73:01: HLA-B*78:01:01: HLA-B*81:01: HLA-B*82:02:01: HLA-C*01:02:01: HLA-C*01:02:11: HLA-C*01:02:29: HLA-C*01:02:30: HLA-C*01:03: HLA-C*01:06: HLA-C*01:08: HLA-C*01:14: HLA-C*01:21: HLA-C*01:30: HLA-C*01:40: HLA-C*02:02:02:01: HLA-C*02:02:02:02: HLA-C*02:10: HLA-C*02:11: HLA-C*02:16:02: HLA-C*02:69: HLA-C*02:85: HLA-C*02:86: HLA-C*02:87: HLA-C*03:02:01: HLA-C*03:02:02:01: HLA-C*03:02:02:02: HLA-C*03:02:02:03: HLA-C*03:03:01: HLA-C*03:04:01:01: HLA-C*03:04:01:02: HLA-C*03:04:02: HLA-C*03:04:04: HLA-C*03:05: HLA-C*03:06: HLA-C*03:100: HLA-C*03:13:01: HLA-C*03:20N: HLA-C*03:219: HLA-C*03:261: HLA-C*03:40:01: HLA-C*03:41:02: HLA-C*03:46: HLA-C*03:61: HLA-C*04:01:01:01: HLA-C*04:01:01:02: HLA-C*04:01:01:03: HLA-C*04:01:01:04: HLA-C*04:01:01:05: HLA-C*04:01:62: HLA-C*04:03:01: HLA-C*04:06: HLA-C*04:09N: HLA-C*04:128: HLA-C*04:161: HLA-C*04:177: HLA-C*04:70: HLA-C*04:71: HLA-C*05:01:01:01: HLA-C*05:01:01:02: HLA-C*05:08: HLA-C*05:09:01: HLA-C*05:93: HLA-C*06:02:01:01: HLA-C*06:02:01:02: HLA-C*06:02:01:03: HLA-C*06:23: HLA-C*06:24: HLA-C*06:46N: HLA-C*07:01:01:01: HLA-C*07:01:01:02: HLA-C*07:01:02: HLA-C*07:01:19: HLA-C*07:01:27: HLA-C*07:01:45: HLA-C*07:02:01:01: HLA-C*07:02:01:02: HLA-C*07:02:01:03: HLA-C*07:02:01:04: HLA-C*07:02:01:05: HLA-C*07:02:05: HLA-C*07:02:06: HLA-C*07:02:64: HLA-C*07:04:01: HLA-C*07:04:02: HLA-C*07:06: HLA-C*07:149: HLA-C*07:18: HLA-C*07:19: HLA-C*07:26: HLA-C*07:30: HLA-C*07:32N: HLA-C*07:384: HLA-C*07:385: HLA-C*07:386: HLA-C*07:391: HLA-C*07:392: HLA-C*07:49: HLA-C*07:56:02: HLA-C*07:66: HLA-C*07:67: HLA-C*08:01:01: HLA-C*08:01:03: HLA-C*08:02:01:01: HLA-C*08:02:01:02: HLA-C*08:03:01: HLA-C*08:04:01: HLA-C*08:112: HLA-C*08:20: HLA-C*08:21: HLA-C*08:22: HLA-C*08:24: HLA-C*08:27: HLA-C*08:36N: HLA-C*08:40: HLA-C*08:41: HLA-C*08:62: HLA-C*12:02:02: HLA-C*12:03:01:01: HLA-C*12:03:01:02: HLA-C*12:08: HLA-C*12:13: HLA-C*12:19: HLA-C*12:22: HLA-C*12:99: HLA-C*14:02:01: HLA-C*14:03: HLA-C*14:21N: HLA-C*14:23: HLA-C*15:02:01: HLA-C*15:05:01: HLA-C*15:05:02: HLA-C*15:13: HLA-C*15:16: HLA-C*15:17: HLA-C*15:96Q: HLA-C*16:01:01: HLA-C*16:02:01: HLA-C*16:04:01: HLA-C*17:01:01:01: HLA-C*17:01:01:02: HLA-C*17:01:01:03: HLA-C*17:03: HLA-C*18:01: HLA-DQA1*01:01:02: HLA-DQA1*01:02:01:01: HLA-DQA1*01:02:01:02: HLA-DQA1*01:02:01:03: HLA-DQA1*01:02:01:04: HLA-DQA1*01:03:01:01: HLA-DQA1*01:03:01:02: HLA-DQA1*01:04:01:01: HLA-DQA1*01:04:01:02: HLA-DQA1*01:05:01: HLA-DQA1*01:07: HLA-DQA1*01:10: HLA-DQA1*01:11: HLA-DQA1*02:01: HLA-DQA1*03:01:01: HLA-DQA1*03:02: HLA-DQA1*03:03:01: HLA-DQA1*04:01:02:01: HLA-DQA1*04:01:02:02: HLA-DQA1*04:02: HLA-DQA1*05:01:01:01: HLA-DQA1*05:01:01:02: HLA-DQA1*05:03: HLA-DQA1*05:05:01:01: HLA-DQA1*05:05:01:02: HLA-DQA1*05:05:01:03: HLA-DQA1*05:11: HLA-DQA1*06:01:01: HLA-DQB1*02:01:01: HLA-DQB1*02:02:01: HLA-DQB1*03:01:01:01: HLA-DQB1*03:01:01:02: HLA-DQB1*03:01:01:03: HLA-DQB1*03:02:01: HLA-DQB1*03:03:02:01: HLA-DQB1*03:03:02:02: HLA-DQB1*03:03:02:03: HLA-DQB1*03:05:01: HLA-DQB1*05:01:01:01: HLA-DQB1*05:01:01:02: HLA-DQB1*05:03:01:01: HLA-DQB1*05:03:01:02: HLA-DQB1*06:01:01: HLA-DQB1*06:02:01: HLA-DQB1*06:03:01: HLA-DQB1*06:09:01: HLA-DRB1*01:01:01: HLA-DRB1*01:02:01: HLA-DRB1*03:01:01:01: HLA-DRB1*03:01:01:02: HLA-DRB1*04:03:01: HLA-DRB1*07:01:01:01: HLA-DRB1*07:01:01:02: HLA-DRB1*08:03:02: HLA-DRB1*09:21: HLA-DRB1*10:01:01: HLA-DRB1*11:01:01: HLA-DRB1*11:01:02: HLA-DRB1*11:04:01: HLA-DRB1*12:01:01: HLA-DRB1*12:17: HLA-DRB1*13:01:01: HLA-DRB1*13:02:01: HLA-DRB1*14:05:01: HLA-DRB1*14:54:01: HLA-DRB1*15:01:01:01: HLA-DRB1*15:01:01:02: HLA-DRB1*15:01:01:03: HLA-DRB1*15:01:01:04: HLA-DRB1*15:02:01: HLA-DRB1*15:03:01:01: HLA-DRB1*15:03:01:02: \
    HLA-DRB1*16:02:01:| $samtools_bin sort --thread $num_processors -m $samtools_sort_memory_per_thread -O BAM - > $sampleid.extract.bam

OUT=$?
if [ ! $OUT -eq 0 ];then
    echo 'Something went wrong while running bwa/samtools to align extracted reads to 38DH_NoAlt (38DH_NoAlt)'
    exit 1
fi

#rm $sampleid.tmp.extract*

echo ">>>>>>>>>>>>>> indexing extracted bam (38DH_NoAlt)"
$samtools_bin index $sampleid.extract.bam

echo ">>>>>>>>>>>>>> bamUtil fastq extraction (38DH_NoAlt)"
$bamUtil bam2FastQ --in $sampleid.extract.bam --gzip --firstOut $sampleid\_extract_1.fq.gz --secondOut $sampleid\_extract_2.fq.gz --unpairedOut $sampleid\_extract.unpaired.fq.gz &> /dev/null

OUT=$?
if [ ! $OUT -eq 0 ];then
    echo '$bamUtil fastq extraction Failed! (38DH_NoAlt)'
    exit 1
else
    rm $sampleid.extract.bam* $sampleid\_extract.unpaired.fq.gz
fi

echo ">>>>>>>>>>>>>> bwa mem to hla panel for Kourami "
$bwa_bin mem -t $num_processors $merged_hla_panel $sampleid\_extract_1.fq.gz $sampleid\_extract_2.fq.gz | $samtools_bin view -Sb - > $bam_for_kourami
OUT=$?
if [ ! $OUT -eq 0 ];then
    echo 'bwa alignment of extracted reads to HLA panel faild...'
    exit 1
fi
