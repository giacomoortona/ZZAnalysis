#!/bin/bash

if [ "$1" == "" ]; then
 echo "Specify set name, e.g. 140125"
 exit 1;
fi

if [ "$2" == "" ]; then
    DESTSET=$1
else
    DESTSET=$2
fi
prodname="rootuples/${1}/PRODFSR/"
dirname="trees/${DESTSET}/PRODFSR/"

#echo $prodname
#echo $dirname

if [ -e $dirname ]; then
    rm $dirname/data/*
    rm $dirname/4mu/*
    rm $dirname/4e/*
    rm $dirname/2mu2e/*
    rm $dirname/CR/*
else
    mkdir -p $dirname/data
    mkdir $dirname/4mu
    mkdir $dirname/4e
    mkdir $dirname/2mu2e
    mkdir $dirname/CR
fi

while read i; do 
    if [[ "$i" != *"#"* ]] && [[ "$i" != "" ]];
    then
	if [ "$i" == "data" ]; then
	    ./run_HZZ4l 0  0 DoubleMu  $prodname/ZZ4lAnalysis_DoubleMu.root  $dirname/data/HZZ4lTree_DoubleMu.root
	    ./run_HZZ4l 1  0 DoubleEle $prodname/ZZ4lAnalysis_DoubleEle.root $dirname/data/HZZ4lTree_DoubleEle.root
	    ./run_HZZ4l 2  0 DoubleOr  $prodname/ZZ4lAnalysis_DoubleOr.root  $dirname/data/HZZ4lTree_DoubleOr.root
	    ./run_HZZ4l_CR 0 DoubleOr  $prodname/ZZ4lAnalysis_DoubleOr.root  $dirname/CR/HZZ4lTree_DoubleOr_CRZLLTree.root
            #./run_HZZ4l_CR 0 DoubleMu  $prodname/ZZ4lAnalysis_DoubleMu.root  $dirname/CR/HZZ4lTree_DoubleMu_CRZLLTree.root
            #./run_HZZ4l_CR 0 DoubleEle  $prodname/ZZ4lAnalysis_DoubleEle.root $dirname/CR/HZZ4lTree_DoubleEle_CRZLLTree.root
	else
	    FullName=`echo $i | sed s/"ZZ4lAnalysis_"//`
	    ./run_HZZ4l 0  0 $i $prodname/${i}.root $dirname/4mu/HZZ4lTree_$FullName.root;
	    ./run_HZZ4l 1  0 $i $prodname/${i}.root $dirname/4e/HZZ4lTree_$FullName.root;
	    ./run_HZZ4l 2  0 $i $prodname/${i}.root $dirname/2mu2e/HZZ4lTree_$FullName.root;
#	    ./run_HZZ4l_CR 0 $i $prodname/${i}.root $dirname/CR/HZZ4lTree__CRZLLTree.root
	fi
    fi
done < Samples_7TeV.txt
