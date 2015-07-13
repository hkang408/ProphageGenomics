def phage_map(path, genome):    #h = hash of gene names
    counter = 0
    switch = 0
    pp_counter = 0
    try:
        a = open(str(path)+'/prophage_tbl.txt', 'rU')
    except:
        print str(path) + ':file not found'
        return

    irofile = iter(a)
    HASH = {}
    a.readline()
    for line in a:
        i = line.split('\t')
        name = str(i[1]).rstrip()
        name = name.upper()

        if switch == 0: # not in pp
            if i[9].startswith('0'): ## Is not a phage
                pass
            elif i[9].startswith('1'): ## Beginning of phage
                MAP  = []
                switch = 1
                pp_counter +=1
                contig = i[2].rstrip()
                MAP.append(ID(name))
        
        elif switch == 1:
            if (i[9].startswith('1') and contig == i[2]): #Same prophage continued.
                pp_counter += 1
                MAP.append(ID(name))
                
            elif (i[9].startswith('1') and contig != i[2]):   # contigs changed in pp. .: 2 prophages
                counter += 1
                contig = i[2]                          
                if int(pp_counter) < 3:
                    MAP = []
                    counter -= 1
                    pp_counter = 1
                    pass
                elif int(pp_counter) >= 350: ### LOOK TO CHANGE AFTER CONSULTING WITH KATE
                    MAP = []
                    counter -= 1
                    pp_counter = 1
                else:
                    HASH[str(genome)+'_'+str(counter)] = MAP
                    MAP = []
                    MAP.append(ID(name))
                    pp_counter = 1
            elif i[9].startswith('0'):
                counter += 1
                switch = 0
                if int(pp_counter) < 3:
                    MAP = []
                    counter -= 1
                    pp_counter = 1
                elif int(pp_counter) >= 350: ### LOOK TO CHANGE AFTER CONSULTING WITH KATE
                    MAP = []
                    counter -= 1
                    pp_counter = 1
                else:
                    HASH[str(genome)+'_'+str(counter)] = MAP
                    MAP = []
                    pp_counter = 0
    a.close()
    return HASH

def automate():
    a = open('genomesNR.txt','U')
    b = open('genomeNR_heatmap.txt', 'w')

    for i in a:
        i = i.strip()
        i = i.replace('\n','')
        path = '/home3/katelyn/KBase/prophages/SEED/'+str(i)
        
        try:
            HASH = phage_map(path,i)
            for pp in HASH:
                LIST = HASH[pp]     # LIST contains array of each gene name encountered in prophage region
                b.write(str(pp))
                for each in LIST:
                    b.write('\t' + str(each))
                b.write('\n')
                    
        except:
            pass
    a.close()
    b.close()

def ID(name):
    if 'INTEGRASE' in name:
        return 'A'
    elif 'PORTAL' in name:
        return 'B'
    elif 'TAIL' in name:
        return 'C'
    elif 'HOLIN' in name:
        return 'D'
    elif 'CAPSID' in name:
        return 'E'
    elif 'ENDONUCLEAS' in name:
        return 'F'
    elif 'TRANSPOSAS' in name:
        return 'G'
    elif 'BASEPLATE' in name:
        return 'H'
    elif 'PROTEASE' in name:
        return 'I'
    elif 'LYSIN' in name:
        return 'J'
    elif 'TERMINASE' in name:
        return 'K'
    elif 'DNA POLYMERASE' in name:
        return 'L'
    elif 'SINGLE STRANDED DNA BINDING' in name.replace('-',' '):
        return 'M'
    elif 'HYPOTHETICAL' in name:
        return 'N'
    else:
        return 'O'
'''
'''
def pick(array):
    if 'A' in array:
        if 'B' in array:
            if 'C' in array:
                if 'D' in array:
                    if 'E' in array:
                        return True
    else:
        return False

def filter_heatmap():
    a = open('genomeNR_heatmap.txt', 'rU')
    b = open('5_heatmapNR.txt', 'w')
    for line in a:
        line = line.rstrip()
        i = line.split('\t')
        genes = i[1:]
        check = pick(genes)
        if check == True:
            int_ = (genes.index('A')+1.0)/len(genes)
            if int_ < 0.5:
                b.write(line+'\n')
            elif int_ >= 0.5:
                b.write(i[0])
                for each in reversed(genes):
                    b.write('\t'+str(each))
                b.write('\n')
        elif check == False:
            pass

def reorder():
    a = open('5_heatmapNR.txt', 'rU')
    b = open('INTfront_heatmapNR.txt', 'w')
    for line in a:
        i = line.split('\t')
        genes = i[1:]
        if genes.index('A') == 0:
            b.write(line)
        elif genes.index('A') != 0:
            int_pos = genes.index('A')
            left = genes[0:int_pos]
            right = genes[int_pos:]
            new_genes = right + left
            b.write(i[0])
            for each in new_genes:
                each = each.rstrip()
                b.write('\t'+each)
            b.write('\n')
    a.close()
    b.close()

automate()
filter_heatmap()
reorder()

def recircularize(): # Will reorder genes to ensure tail goes at end. 
    a = open('NR5_heatmap.txt', 'rU')
    b = open('ASM_heatmap.txt', 'w')
    H = {}
    for line in a:
        line = line.rstrip()
        i = line.split('\t')
        SEQ = i[1:]
        tail_pos = 0
        cap_pos = 0
        tail_count = 0
        cap_count = 0
        length = len(SEQ)
        ID = i[0].rstrip()
        H[ID] = length
        for idx,gene in enumerate(SEQ):
            if gene == 'C':
                if tail_count == 0:
                    cut_here = idx
                pos = (float(idx)*1.0)/len(SEQ)
                tail_pos = tail_pos + pos
                tail_count +=1
     #           print tail_pos
            elif gene == 'B':
                pos = (float(idx)*1.0)/len(SEQ)
                cap_pos = cap_pos + pos
                cap_count +=1
     #           print 'cap' + str(cap_pos)
    
        tail_pos = tail_pos/tail_count
        cap_pos = cap_pos/cap_count
    
        if tail_pos > cap_pos:
            b.write(line+'\n')
        elif tail_pos < cap_pos:
            right = SEQ[cut_here:]
            left = SEQ[:cut_here]
            left.reverse()
            right.reverse()
            J = left +right
            b.write(ID)
            for k in left:
                b.write('\t'+k)
            for k in right:
                b.write('\t'+k)
            b.write('\n')
            
            
    
