'''
The basical approach of this project is as follows:

    1. Genes annotated by automated software (e.g. RAST) can occasionally incorrectly denote a functional gene
        by % homology despite missing a significat section of the gene due to a phage insertion.
    2. The improper annotations often occur outside of phage regions, so I take 4 genes on either
        side (8 genes total per propaphage) and generate new fasta files for them.
    3. These fasta files are BLASTed againsted the ncbi database to find occurances of high % identity and
        low % coverage. These occurences are suspected to be insertion sites.
    4. The 'original' functional gene hit towards NCBI will be then BLASTED against the two prophage regions
        with the hope that it would identify parts of each gene on either side.
    
'''


def insertion_genes(genome):            # Used to generate the gene IDs of the 4 genes at the edge of each prophage region.  


    pp_counter = 0 # Number of genes in prophage region
    counter = 0 # Prophage number designator
    counter2 = 0 # To grab the two genes ahead
    switch = 0 # pp designator
    switch2 = 0 # found pp designator
    path = '/home3/katelyn/KBase/prophages/SEED/' + str(genome)
    HASH = {}
    flag = 0 # Counter for initializing array
    
    try:
        a = open(path+'/prophage_tbl.txt', 'rU')
    except:
        print 'file not found'
        return
        
    INITIAL = []
    LeftHolder = [] # This is where the left region will be held before passed to hash once right region clears. 
        
    for line in a:
        i = line.split('\t')
        flag += 1
        ID = i[0]
        INITIAL.append(ID)
        if flag > 4:        # Limits INITIAL to 4 elements. 
            INITIAL.pop(0)

        if switch == 0:
            if i[9].startswith('0'): # Is not a phage
                counter2+=1
                if (counter2 == 2 and switch2 == 1): # End of prophage from previous 2 lines
                    counter +=1
                    if (pp_counter > 5 and contig == i[2]): 
                        RightHolder = INITIAL
                        counter2 = 0
                        HASH[str(genome)+'_pp_'+str(counter)+'\t'+str(contig)] = LeftHolder+RightHolder
                        LeftHolder = []
                        RightHolder = []
                    elif (pp_counter < 5 or contig != i[2]):
                        counter -=1
                        LeftHolder = []

                    counter2 = 0
                    switch2 = 0
                    pp_counter = 0
                    
                else:
                    pass

            elif i[9].startswith('1'): # Beginning of phage
                counter2 = 1            # Position counter  (Works for non-phage regions too [See line: 36])
                switch = 1              # In pp region
                pp_counter = 1
                contig = i[2]

        elif switch == 1:
            if (i[9].startswith('1') and contig == i[2]): #Same prophage continued.
                pp_counter +=1
                counter2 +=1
                if pp_counter == 2:
                    LeftHolder = LeftHolder + INITIAL

            elif (i[9].startswith('1') and contig != i[2]):  # Contig cutoff
                counter += 1
                counter2 = 1
                
                if pp_counter > 5:
                    RightHolder = INITIAL[1:-1]
                    HASH[str(genome)+'_pp_'+str(counter)+'\t'+str(contig)] = LeftHolder + RightHolder
                    LeftHolder = []
                    RightHolder = []
                    
                elif pp_counter < 5:
                    counter -= 1
                    contig = i[2]                   
                RightHolder = []
                INITIAL = INITIAL[-1:]
                flag = 1
                pp_counter = 1
                contig = i[2]

            elif (i[9].startswith('0') and contig == i[2]):
                counter2 = 1
                switch = 0
                switch2 = 1

    return HASH


def automate():             # Used to locate gene ID of the genes of interest and write files per genome. 
    a = open('complete_genomes.txt', 'rU')
    for line in a:
        line = line.rstrip()
        b = open(str(line)+'_gene.txt', 'w')
        HASH = insertion_genes(line)
        for key in HASH:
            b.write(key)
            for genes in HASH[key]:
                b.write('\t'+genes)
            b.write('\n')
        b.close()





def parser(ID):                 # This function is used to generate the fasta files that will be BLASTED
    a = open('/home3/hkang408/gene_insertions/'+str(ID), 'rU')


    genome = ID[:ID.find('_gene')]
    b = open(str(genome)+'_gene.fasta', 'w')
    
    for line in a:
        i = line.split('\t')
        contig = i[1].rstrip()
        pp = i[0].rstrip()
        genome = pp[:pp.find('_')]
        number = pp[pp.find('pp_')+3:]
        for each in i[2:]:
            each = each.rstrip()
            seq = get_seq(genome,each,contig)
            b.write('>'+str(each)+'_'+str(number)+'\n'+str(seq)+'\n')

    b.close()
    
def get_seq(_genome_,_ID_,_contig_):            # Retrieves AA sequence of file
    import Bio
    from Bio import SeqIO
    handle = open('/home3/katelyn/KBase/prophages/SEED/'+str(_genome_)+'/'+str(_genome_)+'.gbk')
    parse = SeqIO.parse(handle, 'genbank')
    records = list(parse)
    STR =''

    for contig in records:
        if contig.name == _contig_:
            for i in range(0, len(contig.features)):
                try:
                    ID =  contig.features[i].qualifiers['locus_tag'][0]
                    if str(ID) == str(_ID_):
                        STR = contig.features[i].qualifiers['translation'][0]
                    
                except:
                    pass
    handle.close()
    
    return STR
            

def write_fastas():
    import os
    directory = os.listdir('/home3/hkang408/gene_insertions/')
    for i in directory:
        if i.startswith('g.'):
            parser(i)


def _query(argA, argB):
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    pipe = Popen(['blastn',
    '-query', argA,
    '-subject', argB,
    '-outfmt', '6',
    '-word_size', '7',
    '-dust', 'no'],
    stdin=PIPE,
    stdout=PIPE)
    print pipe.communicate()[0]
