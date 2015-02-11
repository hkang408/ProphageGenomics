'''
Creates fasta files of all integrase genes found in prophage region. 

'''


def grab_integrase(path,name):    #h = hash of gene names
    counter = 0
    switch = 0
    pp_counter = 0
    try:
        a = open(str(path)+'/prophage_tbl.txt', 'rU')
    except:
        print str(path) + ':file not found'
        return

    INT = []
    
    HASH = {}
    for line in a:
        i = line.split('\t')
        combo = (str(i[1]),str(i[0]), str(i[2]))

        if switch == 0: # not in pp
            if i[9].startswith('0'): ## Is not a phage
                pass
            elif i[9].startswith('1'): ## Beginning of phage
                MAP  = []
                switch = 1
                pp_counter +=1
                contig = i[2]
                MAP.append(combo)
        
        elif switch == 1:
            if (i[9].startswith('1') and contig == i[2]): #Same prophage continued.
                pp_counter += 1
                MAP.append(combo)
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

    GENEID = []
    for key in HASH:
        for arr in HASH[key]:
            if 'INTEGRASE' in arr[0].upper():

                GENEID.append(arr[1:])
                

    return GENEID


def automate():
    a = open('complete_genomes.txt','U')

    for i in a:
        b = open('compl_genome_int.txt', 'a+')

        i = i.rstrip()
        path = '/home3/katelyn/KBase/prophages/SEED/'+str(i)
        
        integrases = grab_integrase(path,i)
        
        for each in integrases:
            SEQ = get_seq(i,each[0],each[1])
            b.write('>'+each[0]+'\n'+str(SEQ))
            b.write('\n')
        b.close()                    
    a.close()

def get_seq(_genome_,_ID_,_contig_):
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

