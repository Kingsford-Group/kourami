#!/bin/bash

pushd `dirname $0` > /dev/null
SCRIPTD=`pwd`
popd > /dev/null

jvm_memory=4G
resource_dir=$SCRIPTD/../resources
db_base=$SCRIPTD/../custom_db
#input_msa=$SCRIPTD/IMGT/alignments
imgt_ver_num=0
me=`basename $0`

kourami=$SCRIPTD/../target/Kourami.jar
if [ ! -e "$SCRIPTD/../target/Kourami.jar" ];then
    if [ -e "$SCRIPTD/../build/Kourami.jar" ];then
	kourami=$SCRIPTD/../build/Kourami.jar
    elif [ -e "$SCRIPTD/../Kourami.jar"];then
	echo "Could not find Kourami.jar. The jar must be located under target, build or kourami installation directory."
	exit 1
    fi
fi

function usage {
    echo "IMGT/HLA DB formatter for Kourami"
    echo "Note: Run this script to generate Kourami-formatted IMGT/HLA gene-wise MSA"
    echo "      and the reference panel sequences for read alignment"
    echo
    echo "USAGE: <PATH-TO>/$me -i [IMGT/HLA alignments dir] <optional parameters>"
    echo
    echo "------------------ Required Parameters -----------------"
    echo " -i [input_dir]   : path to IMGT/HLA alignments directory "
    echo " -n [hla_nom_g]   : path to matching hla_nom_g.txt file to "
    echo "                    the input alignments directory."
    echo
    echo "------------------ Optional Parameters -----------------"
    echo " -v [ver_number]  : version number is automatically taken from "
    echo "                    IMGT/HLA alignment files. If specified, it "
    echo "                    overrides the version in IMGT/HLA alignment"
    echo "                    files."
    echo " -o [output_dir]  : name of the directory the output will be "
    echo "                  : written to. Output is written to [ver_number]"
    echo "                    directory under [output_dir]. "
    echo "                    (default : custom_db under Kourami installation.)"
    echo " -h               : print this message."
    echo
    exit 1
}

# print usage when no argument is given
if [ $# -lt 1 ]; then
    usage
fi

while getopts i:n:v:o:h FLAG; do
    case $FLAG in 
	i) 
	    input_msa=$OPTARG
	    ;;
	n)
	    nomg=$OPTARG
	    ;;
	v)
	    if [ "$OPTARG" == "0" ]; then
		echo 'Invalid version number: $OPTARG'
		usage
	    fi
	    imgt_ver_num=$OPTARG
	    ;;
	o)
	    db_base=$OPTARG
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

if [[-z "$input_msa" || -z "$nomg"]];then
    echo "missing required parameters"
    usage
fi


if [ ! -e "$nomg" ];then
    echo "$nomg could NOT be found."
    usage
fi

shift $((OPTIND-1))

if [ $# -gt 0 ]; then
    echo "Unrecognized paramester(s) $@"
    usage
fi

if [ ! -e "$resource_dir/HLA_decoys.fa" ];then
    echo "Missing decoy sequence for Kourami panel in the resrouce directory. Please git pull or clone"
    exit 1
fi

if [ ! -e "$resource_dir/DRB5_gen.txt" ];then
    echo "Missing DRB5_gen.txt in the resource directory. Please git pull or git clone"
    exit 1
else
    cp $resource_dir/DRB5_gen.txt $input_msa/.
fi

mkdir -p $db_base
OUT=$?
if [ ! $OUT -eq 0 ];then
    echo "  Cannot create $db_base."
    exit 1
fi

logfilen=`(date +kourami_formatIMGT.%H%M%S%m%d%y.log)`
logfile=$db_base/$logfilen

echo ">>>>>>>>>>>>>> IMGT/HLA DB --> Kourami formatted DB/panel"
java -Xmx$jvm_memory -cp $SCRIPTD/../target/Kourami.jar FormatIMGT $input_msa $imgt_ver_num $db_base 2> $logfile

OUT=$?

#### if something has gone wrong during formatting
if [ ! $OUT -eq 0 ];then
    echo "    An error has occurred while formating IMGT/HLA DB."
    echo "    See log file: $logfile"
    exit 1
else
    if [ "$imgt_ver_num" == "0" ];then
	echo "Getting IMGT/HLA DB release number automatically"
	verline=`head -n1 $logfile`
	if [[ $verline == IMGTver* ]];then
	    IFS=' ' read -a tokens <<< "${verline}"
	    echo IMGT/HLA Release: ${tokens[1]}
	    imgt_ver_num=${tokens[1]}
	fi
    else
	echo "Using the user-input version: $imgt_ver_num"
    fi
    finalLogFile=$db_base/$imgt_ver_num/$logfilen
    if [ -e $db_base/$imgt_ver_num ];then
	mv $logfile $finalLogFile
	echo "Formatting finished. (logfile: $finalLogFile)"
    else
	echo "Formatting finished. (logfile: $logFile)"
    fi
    echo
    outdir=`readlink -e $db_base/$imgt_ver_num/`
    echo "-------------------------------------------------"
    echo " Kourami Formatted db written to: "
    echo " $outdir"
    echo "-------------------------------------------------"
    echo
fi


#### putting the panel sequence together
touch $db_base/$imgt_ver_num/All_FINAL_with_Decoy.fa.gz
panelseq=`readlink -e $db_base/$imgt_ver_num/All_FINAL_with_Decoy.fa.gz`
cat $db_base/$imgt_ver_num/*.merged.fa $resource_dir/HLA_decoys.fa | gzip > $panelseq
cp $nomg $db_base/$imgt_ver_num/.

#### bwa indexing panel sequence
bwa_bin=`(which bwa)`

if [ ! -x "$bwa_bin" ];then 
    echo "[ERROR]: bwa NOT found! Kourami formatted reference panel MUST be indexed before using it."
    echo "Please install the latest copy of bwa and index the panel by running: "
    echo "bwa index $panelseq"
    echo
    exit 1
else
    echo ">>>>>>>>>>>>>> BWA indexing : Kourami reference panel"
    $bwa_bin index $db_base/$imgt_ver_num/All_FINAL_with_Decoy.fa.gz
    echo "-------------------------------------------------"
    echo " Indexed Kourami panel sequence :"
    echo " $panelseq"
    echo "-------------------------------------------------"
    echo
fi
