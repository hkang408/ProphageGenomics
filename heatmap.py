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
    b = open('genome_heatmapNR.txt', 'w')

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
    elif 'TERMINASE' in name:
        return 'F'
    elif 'ENDONUCLEAS' in name:
        return 'G'
    elif 'TRANSPOSAS' in name:
        return 'H'
    elif 'HOLLIDAY' in name:
        return 'I'
    elif 'BASEPLATE' in name:
        return 'J'
    elif 'PROTEASE' in name:
        return 'K'
    elif 'LYSIN' in name:
        return 'L'
    elif 'TOXIC' in name:
        return 'M'
    elif 'LYSOZYM' in name:
        return 'N'
    else:
        return 'O'


'''
'''
def pick(array):
    GENELIST = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O']
    for each in array:
        if each in GENELIST:
            GENELIST.remove(each)
    
    if 'A' not in GENELIST:
        if len(GENELIST) < 10:
            return 1
        else:
            return 0
    else:
        return 0

def filter_heatmap():
    a = open('genome_heatmapNR.txt', 'rU')
    b = open('NR5_heatmap.txt', 'w')
    for line in a:
        line = line.rstrip()
        i = line.split('\t')
        genes = i[1:]
        check = pick(genes)
        if check == 1:
            int_ = (genes.index('A')+1.0)/len(genes)
            if int_ < 0.5:
                b.write(line+'\n')
            elif int_ >= 0.5:
                b.write(i[0])
                for each in reversed(genes):
                    b.write('\t'+str(each))
                b.write('\n')
        elif check == 0:
            pass
        

