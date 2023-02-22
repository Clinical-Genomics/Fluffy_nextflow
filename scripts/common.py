import glob
import os
import math
import numpy

#read the coverage file, store the coverage for each bin together with the gc of each bin
def coverage_tab(file,data):

    for chromosome in data["chromosomes"]:
        data[chromosome]["coverage"] =[]
        data[chromosome]["quality"] = []

    for line in open(file):
        if not line[0] == "#":
            content=line.strip().split("\t")
            content[1]=int(content[1])
            content[3]=float(content[3])
            
            if not content[0] in data["chromosomes"]:
                continue
                
            if content[0] in data:
                data[content[0]]["coverage"].append(content[3])
                if len(content) > 4:
                    data[content[0]]["quality"].append(float(content[4]))

    for chromosome in data["chromosomes"]:
        data[chromosome]["coverage"]=numpy.array(data[chromosome]["coverage"])
        data[chromosome]["quality"]=numpy.array(data[chromosome]["quality"])

    return data


#read the gc content-tab file
def gc_tab(file):
    data={}
    bin_size=-1
    data["chromosomes"] = []
    
    for line in open(file):
        if not line[0] == "#":
            content=line.strip().split("\t")
            content[3]=float(content[3])
            if not "bin_size" in data:
                data["bin_size"] = int(content[2])-int(content[1])
            if not content[0] in data["chromosomes"]:
                data["chromosomes"].append(content[0])
                
            if not content[0] in data:
                data[content[0]] ={}
                data[content[0]]["GC"] =[]
           
            data[content[0]]["GC"].append(content[3])

    for chromosome in data["chromosomes"]:
        data[chromosome]["GC"]=numpy.array(data[chromosome]["GC"])

    return data

#compute a gc histogram
def gc_hist(data,coverage_cutoff,size_cutoff,refQ):
    gc_dictionary={}
    #mean_values={}
    #std={}
    #for chromosome in data["chromosomes"]:
    #    mean_values[chromosome] = numpy.mean(data[chromosome]["coverage"])
        
    for chromosome in data["chromosomes"]:
        for i in range(0,len(data[chromosome]["coverage"])):
            if "Y" in chromosome or "y" in chromosome or "X" in chromosome or "x" in chromosome:
                pass
            else:
                if not data[chromosome]["GC"][i] in gc_dictionary:
                    gc_dictionary[ data[chromosome]["GC"][i]  ]=[]
                
                if data[chromosome]["coverage"][i] > 0:
                    if len(data[chromosome]["quality"]) != 0:
                        if data[chromosome]["quality"][i] >= refQ:
                            gc_dictionary[ data[chromosome]["GC"][i]  ].append(data[chromosome]["coverage"][i])
                    else:
                        gc_dictionary[ data[chromosome]["GC"][i]  ].append(data[chromosome]["coverage"][i])
    hist={}
    for gc in gc_dictionary:
        hist[gc]=[-1,-1]
        if len(gc_dictionary[gc]) > size_cutoff:
            bin_coverage= numpy.median(gc_dictionary[gc])
            if bin_coverage < coverage_cutoff:
                hist[gc]=[bin_coverage,len(gc_dictionary[gc])]

    if not hist:
        print ("Error: Too many low quality regions! consider rerunning the analysis using a smaller --size_cutoff, and  less strict regions masking")
        quit()

    return(hist)


#retrieve the data from the selected region
def regional_cn_est(Data,GC_hist,region,Q):
    chromosome=region[0]
    start=int(region[1])
    end=int(region[2])
    if not chromosome in Data:
        chromosome=chromosome.replace("chr","")

    CN_list=[]
    GC_list=[]
    REF_list=[]

    bins = 0
    used_bins=0
    
    bin_size=Data["bin_size"]
    
    element = 0
    pos= int(math.floor(start/float(bin_size))*bin_size)
    nextpos = pos + bin_size
    while(nextpos > start and pos < end and pos/bin_size < len(Data[chromosome]["coverage"]) ):

            bins += 1
            element=int(pos/bin_size);
            bin_GC=Data[chromosome]["GC"][element]
            bin_coverage=Data[chromosome]["coverage"][element]
            Pass= True
            if len(Data[chromosome]["quality"]) != 0:
                if Data[chromosome]["quality"][element] < Q:
                    Pass = False
            
            if not bin_GC == -1 and not GC_hist[bin_GC][0] == -1:
                if Pass:
                    used_bins += 1
                    CN_list.append(bin_coverage/GC_hist[bin_GC][0])
                    GC_list.append(bin_GC)
                    REF_list.append(GC_hist[bin_GC][0])

            pos+=bin_size;
            nextpos=pos+bin_size;

    if used_bins > 0:
        CN=numpy.median(CN_list)
        GC=numpy.median(GC_list)
        ref=numpy.median(REF_list)
    else:
        CN=-1
        GC=-1
        ref=-1
        
        ci=-1
        
    return([CN,GC,end-start,ref,bins,used_bins,CN_list])


    
