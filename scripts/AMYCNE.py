import argparse
import glob
import os
import common
import sys

import numpy


parser = argparse.ArgumentParser("AMYCNE a copy number estimation toolkit",add_help=False)
parser.add_argument('--genotype' , action="store_true" ,help="compute the copy number in selected regions")
parser.add_argument('--annotate' , action="store_true" ,help="add copy number estimates to structural variant VCF entries")
parser.add_argument('--call' , action="store_true" ,help="perform CNV calling")
parser.add_argument('--hist' , action="store_true" ,help="compute the coverage across each chromosome, return a tab file describing the average coverage, as well as average coverage per contig")
parser.add_argument('--count' , action="store_true" ,help="estimate the copy number of each chromosome")
parser.add_argument('--ff' , action="store_true" ,help="Predict fetal fraction based on the number of reads of chromosome Y")
parser.add_argument('--filt' , action="store_true" ,help="filters the input coverage tab file, prints the filtered and gc corrected version to stdout")
args, unknown = parser.parse_known_args()


if args.ff:

    parser = argparse.ArgumentParser("""AMYCNE-ff: Estimate fetal fraction based on read depth across chromosome Y""")
    parser.add_argument('--ff' , action="store_true" ,help="Predict fetal fraction based on the number of reads of chromosome Y")
    parser.add_argument('--coverage' , type=str, help="the tab file containing coverage tab files")
    parser.add_argument('--Q' , type=int,default=15,help="Minimum average mapping quality of the bins used for copy number estimation default = 15")
    parser.add_argument('--scaling' , type=float,default=2,help="The slope of FFY (default=2)")
    parser.add_argument('--intercept' , type=float,default=0,help="Set the intercept of FFY (default=0)")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--refQ' , type=int,default=10,help="Minimum average mapping quality of the bins used for constructing the reference = 10")
    parser.add_argument('--c_cutoff' , type=int,default=100,help="bins having coverage higher than the cut off value are excluded from the ref calculations")
    parser.add_argument('--s_cutoff' , type=int,default=50,help="bins that have less than the s_cutoff value similar bins are discarded from copy number esitmation")
    args = parser.parse_args()

    Data= common.gc_tab(args.gc)
    Data=common.coverage_tab(args.coverage,Data)
    GC_hist=common.gc_hist(Data,args.c_cutoff,args.s_cutoff,args.refQ)
    autosomal_bins=[]
    y_bins=[]
    x_bins=[]
 
    for chromosome in Data["chromosomes"]:
        if "21" in chromosome or "13" in chromosome or "9" in chromosome or "18" in chromosome or len(chromosome) > 5:
           continue

        for i in range(0,len(Data[chromosome]["coverage"])):
           if not Data[chromosome]["GC"][i] in GC_hist or Data[chromosome]["GC"][i] == -1 or GC_hist[Data[chromosome]["GC"][i]][0] == -1:
               continue

           if Data[chromosome]["coverage"][i] > 0 and Data[chromosome]["quality"][i] < args.Q:  
               continue

           if chromosome == "Y" or chromosome == "chrY": 
               if Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0] > 0.5:
                   continue
               y_bins.append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])
           elif chromosome == "X" or chromosome == "chrX":
               x_bins.append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])

           else:
               autosomal_bins.append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])

    sex="female"
    ff="NA"

    if len(y_bins):
        medY=numpy.median(y_bins)
    else:
        medY=0
    if len(x_bins):
        medX=numpy.median(x_bins)
    else:
        medX=0

    medA=numpy.median(autosomal_bins)

    if medY > 0.005:
        sex= "male"

    print ("sample sex medAutosomal FFY FFX")
    if sex == "male":
        print ("{} {} {} {} {}".format(args.coverage.split("/")[-1].split(".")[0],sex, medA, args.scaling*medY+args.intercept, args.scaling*(1-medX) ) )
    else:
        print ("{} {} {} {} {}".format(args.coverage.split("/")[-1].split(".")[0],sex, medA, args.scaling*medY, args.scaling*(1-medX) ) )
