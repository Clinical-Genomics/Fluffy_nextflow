import argparse
import gzip
               
parser = argparse.ArgumentParser("""Prints tab separated GC content in windows of user set size""")
parser.add_argument('--fa' , type=str, required=True,help="the fasta file")
parser.add_argument('--size' , type=int, default=1000, help="the size of the bins")
parser.add_argument('--replace_chr' , action="store_true", help="remove chr prefix(i.e chr1 -> 1)")
parser.add_argument('--n_mask', action="store_true" , help="regions having more than 10%% content set to N or n is set to -1")
args = parser.parse_args()



def calculate_complexity(section,k):
    complexity=0;
    movement=len(section)-k+1
    pos=0
    
    kmer_combinations={};
    while pos <= movement:
        kmer=section[pos:pos+k]
        if not "N" in kmer:
            if not kmer in kmer_combinations:
                kmer_combinations[kmer]=0
            
            kmer_combinations[kmer] += 1
        pos += 1
        
    total_kmers=0;
    max_kmer=-1;
    for kmer in kmer_combinations:
        total_kmers += kmer_combinations[kmer]
        if kmer_combinations[kmer] > max_kmer:
            max_kmer=kmer_combinations[kmer]
    if total_kmers > 0:
        complexity=max_kmer/float(total_kmers)
    else:
        complexity = -1
        
    if complexity > 0.2:
        complexity == -1    
    return(complexity)


#read the fasta file
reference={}
chromosomes=[]
if not args.fa.endswith(".gz"):
    with open(args.fa, 'r') as f:
        sequence=f.read()
else:
    with gzip.open(args.fa, 'r') as f:
        sequence=f.read()

split_reference=sequence.split(">")
del sequence
del split_reference[0]
#store the reference as a dictionary
for chromosome in split_reference:
    content=chromosome.split("\n",1)
    reference[content[0].strip().split()[0]]=content[1].replace("\n","").upper()
    chromosomes.append( content[0].strip().split()[0])
        
#calculate the coverage of each bin
for chromosome in chromosomes:
    start=0
    end=args.size
    while start < len(reference[chromosome]):
        printed_chromosome=chromosome
        if args.replace_chr:
            printed_chromosome=chromosome.replace("chr","")
        region=reference[chromosome][start:end]
        G=region.count("C")+region.count("G")
        A=region.count("A")+region.count("T")
        N=region.count("N")
        if A+G > 0:
            GC=float(G)/float(A+G)
        elif args.n_mask and N > 0.1:
            GC=-1
        else:
            GC=0
            
        complexity = calculate_complexity(region,3)
        if complexity == -1:
            GC=-1
            
        print("{}\t{}\t{}\t{}".format(printed_chromosome,start,end, round(GC,2) ))
        start=end
        if end+args.size < len(reference[chromosome]):
            end += args.size
        else:
            end=len(reference[chromosome])
