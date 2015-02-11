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
                    

                print 'contig cutoff'
                
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
                if (int(position) > 1000 and len(line) > int(position)+1000):
                    sequence = line[int(position)-1000:int(position)+1000]
                elif (int(position) < 1000 and len(line) > int(position)+1000):
                    sequence = line[0:int(position)+1000]
                elif (int(position) > 1000 and len(line) < int(position)+1000):
                    sequence = line[int(position)-1000:len(line)]
                else:
                    print 'Position setting error'
                    
                switch = 'False'
            else:
                pass
        return str(sequence)
    except:
        return 'nothing'




def auto():
    a = open('genomes.txt', 'rU')
    for line in a:
        line = line.replace('\n','')
        pp_fastas(line)

    a.close()

auto()
        

