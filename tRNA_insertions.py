# Uses the data from the pp_summary file to generate fasta files
# of up to 2000 nt of prophage region flanks. Names each sequence based on genome_prophage#_L/R.

# EST: Aug 14, 2014

def FIND_SEQ(genome, pp_count, contig, start, end):
    import os
    genome = genome.replace('*','')
    genome = genome.replace('\n','')
    genome = genome.strip()
    try:
        path = '/home3/katelyn/KBase/SEED/' + str(genome) + '/contigs'
        a = open(path, 'U')
    except:
        print 'File Missing: ' + str(genome)
    b = open(str(genome) + '_tRNA.fasta', 'a+')
    switch = 0

    for i in a:
        i = i.replace('\n','')
        if i.startswith('>'):
            if str(i) == ('>'+str(contig)):
                switch = 1
            else:
                pass
        else:
            
            if switch == 0:
                pass
            elif switch == 1:
    
                if int(start) < 1000:
                    SEQ_L = i[0:int(start)+1000]
                elif int(start) > 1000:
                    SEQ_L = i[int(start)-1000: int(start)+1000]

                if len(i) < int(end)+1000:
                    SEQ_R = i[int(end)-1000:int(end)]

                elif len(i) > int(end)+1000:
                    SEQ_R = i[int(end)-1000:int(end)+1000]
                switch = 0
    b.write('>' + str(genome) +'_' +str(pp_count)+ '_L\n' + str(SEQ_L) +'\n'
            +'>' + str(genome) + '_' + str(pp_count) + '_R\n' + str(SEQ_R) + '\n')
                
    a.close()
    b.close()

def automate_FIND_SEQ():
    
    a = open('pp_summary.txt', 'U')
    a.readline()
    for line in a:
        i = line.split('\t')
        genome = i[0]
        if '*' in genome:
            genome = genome.replace('*','')
        try:
            pp_count = i[1]
            contig = i[2]
            start = i[3]
            end = i[4]
            FIND_SEQ(genome,pp_count,contig,start,end)
        except:
            print 'ERROR: ' + str(genome)
    a.close()
    b.close()
    

# Identifies tRNAs with identical codons that are present in both flanking regions
# of a genome. 
# Infile:
'''
g.0_1_L 	1	908 	981 	Arg	TCT	0	0	73.38
g.1004_1_L 	1	672 	744 	Phe	GAA	0	0	75.98
g.1004_4_R 	1	1126	1053	Pro	GGG	0	0	47.98
g.1007_1_R 	1	1209	1136	Arg	TCT	0	0	81.13
g.1010_4_L 	1	770 	842 	Lys	CTT	0	0	86.77
g.1010_5_R 	1	1257	1167	SeC(p)	TCA	0	0	71.03

'''
# Outfile:

'''
g.1701_3_R	GCT	119	32
g.1701_3_L	GCT	424	337	
g.23200_4_R	CAT	1545	1617
g.23200_4_L	CAT	1773	1845	
g.23538_1_R	GGA	304	218
g.23538_1_L	GGA	422	336	
g.23702_3_R	GTT	25	98
g.23702_3_L	GTT	734	807

'''
# EST: Aug 30th, 2014
def find_common_tRNA():
    a = open('tRNA.txt', 'rU')
    b = open('tRNAs_present_in_both_flanks.txt', 'w')

    data_L = {}
    data_R = {}
    counter = 0

    for line in a:
        i = line.split('\t')
        genome = i[0].replace(' ','')
        start = i[2].replace(' ','')
        end = i[3].replace(' ','')
        codon = i[5].replace(' ','')
        base = genome[:-1]

        if 'L' in genome:
            if genome not in data_L:
                data_L[genome] = {(start,end):codon}
            elif genome in data_L:
                data_L[genome][(start,end)] = codon
            elif (start,end) in data_R[genome]:
                print 'fuuuuk'

            if base+'R' not in data_R:
                pass
            elif base+'R' in data_R:
                print 'yes'
                for KEY in data_R[base+'R']:
                    CODON = data_L[base+'R'][KEY]
                    if CODON == codon:
                        b.write(str(genome)+'\t' + str(codon) + '\t' + str(start) + '\t' + str(end) + '\n')
                        b.write(str(base+'R')+'\t' + str(CODON)+ '\t')
                        for each in KEY:
                            b.write(each+'\t')
                        b.write('\n')

        elif 'R' in genome: # If the line contains a right flank
            if genome not in data_R:
                data_R[genome] = {(start, end):codon}
            elif genome in data_R:
                data_R[genome][(start,end)] = codon
            elif (start,end) in data_R[genome]:
                print 'fuckkk'

            if base+'L' not in data_L:
                pass
            elif base+'L' in data_L:
                for KEY in data_L[base+'L']:
                    CODON = data_L[base+'L'][KEY]
                    if CODON == codon:
                        b.write(str(genome)+'\t'+str(codon)+'\t'+str(start)+'\t'+str(end)+'\n')
                        b.write(str(base+'L')+'\t'+str(CODON)+'\t')
                        for each in KEY:
                            b.write(each+'\t')
                        b.write('\n')

    a.close()
    b.close()
    return

### Obtains the tRNA sequence from the prophage flank fasta files.
### Infile:
'''
g.1701_3_R	GCT	119	32
g.1701_3_L	GCT	424	337	
g.23200_4_R	CAT	1545	1617
g.23200_4_L	CAT	1773	1845	
g.23538_1_R	GGA	304	218
g.23538_1_L	GGA	422	336	
g.23702_3_R	GTT	25	98
g.23702_3_L	GTT	734	807	
g.23843_1_R	CTT	46	118
g.23843_1_L	CTT	267	339
'''
# To generate file that looks like:
'''
>g.1701_3_R
GGAGAAAGAGGGATTTGAACCCTCGCGCCGGTCACCCGACCTATACCCTTAGCAGGGGCACCTCTTCAGCCTCTTGAGTATTTCTCC
>g.1701_3_L
GGAGAAAGAGGGATTTGAACCCTCGCGCCGGTCACCCGACCTATACCCTTAGCAGGGGCACCTCTTCAGCCTCTTGAGTATTTCTCC
>g.23200_4_R
GCCCTTTAGCTCAGTGGTGAGAGCGAGCGACTCATAATCGCCAGGTCGCTGGTTCAAATCCAGCAAGGGCCA
>g.23200_4_L
GCCCTTTAGCTCAGTGGTGAGAGCGAGCGACTCATAATCGCCAGGTCGCTGGTTCAAATCCAGCAAGGGCCA
>g.23538_1_R
GGAGAGAGAGGGATTCGAACCCTCGATACGGTTTCCCGTATACACGCGTTCCAGGCGTGCGCCTTCAACCACTCGGCCACCTCTCC
>g.23538_1_L
GGAGAGAGAGGGATTCGAACCCTCGATACGGTTTCCCGTATACACGCGTTCCAGGCGTGCGCCTTCAACCACTCGGCCACCTCTCC
>g.23702_3_R
CCGCCATAGCTCAGTTGGTAGTAGCGCATGACTGTTAATCATGATGTCGTAGGTTCGAGTCCTACTGGCGGAG

'''
# EST: Aug 24th, 2014

def get_tRNA_seq(inputseq, start, end):
    file_ = inputseq[:inputseq.find('_')]
    file_ = file_+str('_tRNA.fasta')
    a = open('/home3/hkang408/tRNA/'+str(file_),'rU')

    switch = 0
    for line in a:
        if line.startswith('>'+str(inputseq)):
            switch = 1
            pass
        elif switch == 1:
            fasta = line[int(start):int(end)]
            switch = 0
        elif switch == 0:
            pass
        else:
            pass
    a.close()
    return fasta

def automate_get_tRNA_seq():
    a = open('tRNA_in_both_flanks.txt', 'rU')
    b = open('tRNA_seq.fasta', 'w')
    for i in a:
        i = i.split('\t')
        if int(i[2]) > int(i[3]):
            start = i[3]
            end = i[2]
        elif int(i[3]) > int(i[2]):
            start = i[2]
            end = i[3]
        fasta = get_tRNA_seq(i[0], start, end)
        b.write('>'+str(i[0])+'\n' + fasta+'\n')
        
    a.close()
    b.close()
    return

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

def automate_tRNA_BLAST():
    import os
    import tempfile
    
    a = open('/home3/hkang408/tRNA/OUT2/tRNA_seq.fasta', 'rU')
    tf = tempfile.NamedTemporaryFile()
    counter = 0
    
    for i in a:
        if counter % 2 == 0:        # Even, or 0
            line1 = i
            line1 = line1.replace('\n','')
            counter +=1
        elif counter % 2 != 0:  # Odd, not 0
            tf = tempfile.NamedTemporaryFile()
            line2 = i
            line2 = line2.replace('\n','')
            counter = 0
            tf.write(line1 + '\n' + line2)
            tf.seek(0)
            query = tf.name
            subject = '/home3/hkang408/tRNA/'+str(line1[1:line1.find('_')])+'_tRNA.fasta'
            _query(query, subject)
            
                
    a.close()

def identify():                 #
    a = open('OUTS', 'rU')
    b = open('tRNAs', 'w')
    for line in a:
        i = line.split('\t')
        try:
            query = i[0]
            subject = i[1]
            evalu = i[10]
            if 'L' in query:
                if subject == query[:-1]+'R':
                    if float(evalu) < 0.1:
                        b.write(line)
        except:
            pass


    

