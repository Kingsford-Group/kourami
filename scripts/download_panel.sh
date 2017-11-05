#
# Part of Kourami HLA typer/assembler
# (c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
# See LICENSE for licensing.
#

#!/bin/bash

function print_usage {
    echo "Usage: $0"
    echo ""
    echo "Note: This script downloads Kourami reference panel sequences."
    echo "      The panel is downloaded in '<kourami_intallation>/db'   "
    exit 1;
}

echo "#---------------------------------------------------------------------"
echo "# Run this script to download Kourami panel reference."
echo "#---------------------------------------------------------------------"
echo

pushd `dirname $0` > /dev/null
SCRIPTD=`pwd`
popd > /dev/null
kourami_home=$SCRIPTD/..

urlkouramipanel="https://github.com/Kingsford-Group/kourami/releases/download/v0.9/kouramiDB_3.24.0.tar.gz"

bwa_bin=`(which bwa)`

if [ ! -x "$bwa_bin" ];then
    echo "Please make sure bwa is installed."
    exit 1
fi

if [ $# -eq 0 ]; then
    echo "-----------------------------------------------"
    echo "| Downloading Kourami panel ...               |"
    echo "-----------------------------------------------"
    echo
    wget -O- $urlkouramipanel | tar xzf - -C $kourami_home
    OUT=$?
    if [ ! $OUT -eq 0 ];then
	echo 
	echo "-----------------------------------------------"
	echo "| Could NOT download Kourami reference panel! |"
	echo "-----------------------------------------------"
	exit 1
    fi
    echo "-----------------------------------------------"
    echo "| Indexing the Kourami panel ...              |"
    echo "-----------------------------------------------"
    echo
    $bwa_bin index $kourami_home/db/All_FINAL_with_Decoy.fa.gz 
    OUT=$?
    if [ ! $OUT -eq 0 ];then
	echo 
	echo "-----------------------------------------------"
	echo "| There was a problem indexing Kourami panel  |"
	echo "-----------------------------------------------"
	exit 1
    fi
    echo 
    echo "----------------------------------"
    echo "| Kourami panel is installed at: |"
    echo "----------------------------------"
    readlink_bin=`(which readlink)`
    if [ ! -x "$readlink_bin" ];then
	echo "$kourami_home/db"
    else
	$readlink_bin -f $kourami_home/db
    fi
else
    print_usage
fi

