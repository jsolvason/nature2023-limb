from Bio import SeqIO 
import matplotlib.pyplot as plt
from pyliftover import LiftOver

def faLoadGenome(file_genome):
    '''Loads arbitrary genome'''
    Handle = open(file_genome, "r") 
    chr2seq = SeqIO.to_dict(SeqIO.parse(Handle, "fasta")) 
    chr2seq = {chrom:str(chr2seq[chrom].seq).upper() for chrom in chr2seq}
    return chr2seq

def generate_Chrom2NumStartEnd(chr2seq):
    '''Generates obj for choose_random_chrom_pos()'''
    genomeLen=sum([len(seq) for seq in chr2seq.values()])
    Chrom2NumStartEnd={}  
    n=0
    for chrom,seq in chr2seq.items():
        begin=n
        end=begin+len(seq)
        Chrom2NumStartEnd[chrom]=(begin,end)
        n=end
    return genomeLen,Chrom2NumStartEnd

def choose_random_chrom_pos(genomeLen,Chrom2NumStartEnd):
    '''Chooses random position in the genome.'''
    n=np.random.randint(genomeLen)
    for chrom,(start,end) in Chrom2NumStartEnd.items():
        if n>= start and n<=end:
            return chrom,n-start
        
def write_row(rowList,delim='\t'):
    '''Write a single row of a tsv file.'''
    return delim.join([str(i) for i in rowList])+'\n'

def IupacToAllPossibleSequences(dna):
	'''Takes DNA with IUPAC letters and returns all possible DNA strings with only A/G/T/C.'''
	Round2Seqs={}
	Round2Seqs[0]=[]
	for nt in Iupac2AllNt[dna[0]]:
		Round2Seqs[0].append(nt)

	for i,iupac in enumerate(dna):
		if i==0: continue
			
		Round2Seqs[i]=[]
		
		for seq in Round2Seqs[i-1]:
			for nt in Iupac2AllNt[iupac]:
				thisSeq=seq+nt
				Round2Seqs[i].append(thisSeq)
				
	lastRound=max(Round2Seqs.keys())
	return Round2Seqs[lastRound]
        
def get_kmers(string,k):
	'''Takes DNA sequence as input and a kmer length and YIELDS all kmers of length K in the sequence.'''
	for i in range(len(string)):
		kmer8=string[i:i+k]
		if len(kmer8)==k:
			yield kmer8
            
def revcomp(dna): 
	'''Takes DNA sequence as input and returns reverse complement'''
	inv={'A':'T','T':'A','G':'C','C':'G', 'N':'N','W':'W'}
	revcomp_dna=[]
	for nt in dna:
		revcomp_dna.append(inv[nt])
	return ''.join(revcomp_dna[::-1])

def loadAff(ref):
        '''Loads a preprocessed pbm dataset.'''
        Seq2EtsAff  = {line.split('\t')[0]:float(line.split('\t')[1]) for line in open(ref,'r').readlines()}
        return Seq2EtsAff
    
def mkdir_if_dir_not_exists(out_dir):
    '''Make a directory only if that directory doesnt exist'''
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    
def dprint(d,n=0):
    '''Prints a dictionary'''
    for i,(k,v) in enumerate(d.items()):
        if i<=n:
            print(k,v)
            
def zipdf(df,cols):
    '''iterates through multiple columns of a dataframe'''
    return zip(*[df[c] for c in cols])

def quickfig(x=5,y=5,dpi=150):
    '''quickly generate a matplotlib figure'''
    return plt.subplots(1,figsize=(x,y),dpi=dpi)

# color blind friendly colors
cb={}
cb['lightblue']= [i/255 for  i in [86,180,233]]
cb['green']    = [i/255 for  i in [0,158,115]]
cb['red']      = [i/255 for  i in [213,94,0]]
cb['yellow']   = [i/255 for  i in [240,228,66]]
cb['orange']   = [i/255 for  i in [230,159,0]]
cb['blue']     = [i/255 for  i in [0,114,178]]
cb['pink']     = [i/255 for  i in [204,121,167]]
cb['black']    = [i/255 for  i in [0,0,0]]

# iupac dna definition
Iupac2AllNt= {
        'A':['A'],
        'C':['C'],
        'G':['G'],
        'T':['T'],
        'R':['A','G'],
        'Y':['C','T'],
        'S':['G','C'],
        'W':['A','T'],
        'K':['G','T'],
        'M':['A','C'],
        'B':['C','G','T'],
        'D':['A','G','T'],
        'H':['A','C','T'],
        'V':['A','C','G'],
        'N':['A','C','G','T'],
}

def getLiftoverObject(startGenome,endGenome):
    '''A method to create a liftover object for jsg.liftover()'''
    return LiftOver(startGenome, endGenome)

def read_tsv(fn,pc,header,breakBool=False,sep='\t',pc_list=False):
    '''Read a tsv file'''
    with open(fn,'r') as f:
        
        # If printing columns, skip header
        if header:
            if (pc==True): pass
            else:          next(f)
                
        for i,line in enumerate(f):
            a=line.strip().split(sep)
            if pc:
                if pc_list==False:
                    if i==0:
                        for i,c in enumerate(a):
                            print(i,c)
                        print()
                        if breakBool: break
                        continue
                else:
                    if i==0:
                        print(', '.join([i.replace('-','_').replace(' ','_') for i in a]))
                        if breakBool: break
                        continue
            yield a
            
def flatten_list(l):
    '''Returns [1,2,3,4] for inputted [[1,2],[3,4]]'''
    return [i for sublist in l for i in sublist]

def liftover_pos(liftOverObject,c,s):
    '''A method to liftover a single chrom/pos pair'''
    strloResults = liftOverObject.convert_coordinate(c,s)
    if (strloResults==[]): return False,'NotFound'
    if (strloResults==None): return False,'NotFound'
    else:
        hg38chrom,hg38pos,strand,score=strloResults[0]
        if c!=hg38chrom:     
            return False,'ChrOld≠ChrNew'
        else:                
            return hg38chrom,hg38pos