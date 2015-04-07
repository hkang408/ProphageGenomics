'''
This script will generate nucleotide sequence fasta files of the region surrounding a designated prohage site.
Can be used in conjuncture with both tRNA insertion and gene insertion approaches. 
'''
def pp_fastas(genome):                          
# Setting all necessary parameters
    pp_counter = 0
    counter = 0
    switch = 0 # pp designator
    path = '/home3/katelyn/KBase/prophages/SEED/' + str(genome)
    try:
        a = open(path+'/prophage_tbl.txt', 'rU')
        b = open(str(genome)+'_flanks.fasta', 'a+')
    except:
        print 'file not found'
        return
    irofile = iter(a)

    for line in irofile:
        i = line.split('\t')
        if switch == 0: # not in pp
            if i[9].startswith('0'): ## Is not a phage
                pass
            elif i[9].startswith('1'): ## Beginning of phage
                switch = 1
                pp_counter += 1
                contig = i[2]
                counter += 1
                if int(i[3]) < int(i[4]):               # Designate the 'left' side of pp
                    L_start = i[3]
                elif int(i[3]) > int(i[4]):
                    L_start = i[4]

        elif switch == 1:
            if (i[9].startswith('1') and contig == i[2]): #Same prophage continued.
                pp_counter +=1
                if int(i[3]) < int(i[4]):
                    R_end = i[4]
                elif int(i[3]) > int(i[4]):
                    R_end = i[3]
        
            elif (i[9].startswith('1') and contig != i[2]):   # contigs changed in pp. .: 2 prophages
                if pp_counter > 5:  # Long enough to be considered a prophage
                    sequence1 = RETR_SEQ(path,contig,L_start)
                    sequence2 = RETR_SEQ(path,contig,R_end)
                    b.write('>'+str(genome)+'_'+str(counter)+'|L\n'+str(sequence1)+'\n')
                    b.write('>'+str(genome)+'_'+str(counter)+'|R\n'+str(sequence2)+'\n')
                    contig = i[2]
                    if int(i[3]) < int(i[4]):               # Designate the 'left' side of pp
                        L_start = i[3]
                    elif int(i[3]) > int(i[4]):
                        L_start = i[4]
                    counter += 1
                else:
                    pp_counter = 1
                    counter -= 1
                pp_counter = 1
                
            elif i[9].startswith('0'): # Reached the end of pp
                switch = 0
                if L_start != 0: # as in there was no contig cutoff
                    if pp_counter > 5:  # Long enough to be considered a prophage | at this point all 
                        sequence1 = RETR_SEQ(path,contig,L_start)
                        sequence2 = RETR_SEQ(path,contig,R_end)
                        b.write('>'+str(genome)+'_'+str(counter)+'|L\n'+str(sequence1)+'\n')
                        b.write('>'+str(genome)+'_'+str(counter)+'|R\n'+str(sequence2)+'\n')
                    else:
                        counter -= 1
                pp_counter = 0
                    
    a.close()
    b.close()
    return                   

def RETR_SEQ(path, contig, position):
    try:
        a = open(str(path)+'/contigs', 'rU')
        switch = 'False'
        for line in a:
            line = line.strip()
            line = line.replace('\n','')
            if contig == line[1:]:
                switch = 'True'
            elif switch == 'True':
                if (int(position) > 5000 and len(line) > int(position)+5000):
                    sequence = line[int(position)-6000:int(position)+5000]
                elif (int(position) < 5000 and len(line) > int(position)+5000):
                    sequence = line[0:int(position)+5000]
                elif (int(position) > 5000 and len(line) < int(position)+5000):
                    sequence = line[int(position)-5000:len(line)]
                else:
                    print 'Position setting error'
                    
                switch = 'False'
            else:
                pass
        return str(sequence)
    except:
        return 'nothing'

def auto():
    a = open('genomesNR.txt', 'rU')
    for line in a:
        line = line.replace('\n','')
        pp_fastas(line)
    a.close()

auto()
#####       tRNAscan-SE scripts
def _query(argA):
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    pipe = Popen(['/home3/katelyn/opt/tRNAscan-SE/tRNAscan-SE', '-P',  argA
   ],             #        '-dust', 'no'
    stdin=PIPE,
    stdout=PIPE)
    print pipe.communicate()[0]
    
def tRNAscan():
    a = open('pp_flank_lists.txt','rU')
    path = '/home3/hkang408/pp_flanks_10k_NR/'
    for line in a:
        line = line.rstrip()
        argA = path+str(line)
        _query(argA)
    a.close()




####        Find tRNA fragments on opposite side of prophage
def automate():
    a = open('out_tRNA.txt', 'rU')
    for line in a:
        i = line.split('\t')
        flank = i[0].rstrip()
        tRNAtype = i[4].rstrip()
        start_ = i[2].rstrip()
        end_ = i[3].rstrip()
        if line.startswith('g.'):
            genome = flank[:flank.find('_')]
            if int(start_) < int(end_):
                start = start_
                end = end_
            elif int(start_) > int(end_):
                start = end_
                end = start_
            seq = get_tRNA_seq(genome,flank,start,end)
 #           print seq
            b = open('tmp.txt', 'w')
            b.write('>'+str(flank)+'_'+str(tRNAtype)+'\n'+str(seq))
            b.close()
            _query('tmp.txt', '/home3/hkang408/pp_flanks_10k_NR/'+str(genome)+'_flanks.fasta')           

def get_tRNA_seq(genome,flankID,start,end):
    a = open('/home3/hkang408/pp_flanks_10k_NR/'+str(genome)+'_flanks.fasta', 'rU')
    irofile = iter(a)

    for line in a:
        if flankID in line:
            fullseq = next(irofile)
            seq = fullseq[int(start):int(end)]
    return seq

def _query(argA, argB):                 # To blast for tRNA fragments
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

#automate()

def analyze_blast_out():
    a = open('OUT', 'rU')
    for line in a:
        try:
            i = line.split('\t')
            evalue = i[10]
            query = i[0]
            subject = i[1]
            base = query[:query.find('|')]
            if '|L' in query:
                if subject == base+'|R':
                    if float(evalue) < 0.01:
                        print line.rstrip()
            elif '|R' in query:
                if subject == base+'|L':
                    if float(evalue) < 0.01:
                        print line.rstrip()
        except:
            pass
    a.close()

