def word_count(gene, h):
    i = gene.split(' ')
    for word in i:
        if word not in h:
            h[word] = 1
        elif word in h:
            value = int(h[word])
            new_val = value + 1
            h[word] = new_val
    return h

def phage_map(path,name):    #h = hash of gene names
    counter = 0
    switch = 0
    pp_counter = 0
    try:
        a = open(str(path)+'/prophage_tbl.txt', 'rU')
    except:
 #       print str(path) + ':file not found'
        return

    irofile = iter(a)
    HASH = {}
    for line in irofile:
        i = line.split('\t')
        if switch == 0: # not in pp
            if i[9].startswith('0'): ## Is not a phage
                pass
            elif i[9].startswith('1'): ## Beginning of phage
                MAP  = []
                switch = 1
                pp_counter +=1
                contig = i[2]
                MAP.append(str(i[1]))
        
        elif switch == 1:
            if (i[9].startswith('1') and contig == i[2]): #Same prophage continued.
                pp_counter += 1
                MAP.append(str(i[1]))
                contig = i[2]
            elif (i[9].startswith('1') and contig != i[2]):   # contigs changed in pp. .: 2 prophages
                counter += 1
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
                    HASH[str(name)+'_'+str(counter)] = MAP
                    MAP = []
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
                    HASH[str(name)+'_'+str(counter)] = MAP
                    MAP = []
                    pp_counter = 1
    a.close()
    return HASH

def automate():
    a = open('whole_genomes.txt','U')
    b = open('whole_genome_pp.txt', 'w')

    for i in a:
        i = i.strip()
        i = i.replace('\n','')
        path = '/home3/katelyn/KBase/prophages/SEED/'+str(i)
        
 #       try:
        HASH = phage_map(path,i)
        for pp in HASH:
            LIST = HASH[pp]     # LIST contains array of each gene name encountered in prophage region
            ORDER = genemap(LIST)
            b.write(str(pp)+'\t')
            for j in ORDER:
                b.write(str(j)+'\t')
            b.write('\n')
                    
#        except:
 #           pass
    a.close()
    b.close()

def genemap(LIST):
    gene_dictionary = {
    'INTEGRASE':1,
    'PORTAL':2,
    'TAIL': 3,
    'HOLIN':4,
    'TRANSPOSA':5,
    'TERMINA':6,
    'ENDONUCLEAS':7,
    'CAPSID':8,
    'HOLLIDAY':9,
    'BASEPLATE':10,
    'PROTEASE':11,
    'LYSIN':12,
    'TOXI':13,
    'LYSOZYM':14}
    SEQ = []
    for gene in LIST:
        flag = 0
        counter = 1
        for key in gene_dictionary.keys():
            if flag == 1: # Key had hit
                pass
            elif flag == 0:
                if key in gene.upper():
                    SEQ.append(gene_dictionary[key])
                    flag = 1
                else:
                    counter +=1
            if counter == 15:
                SEQ.append(0)
                counter +=1
    return SEQ

