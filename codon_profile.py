def codon_profile(file_):
    a = open('/home3/hkang408/pp_genes/'+str(file_))
    genome = file_[:file_.find('.fasta')]
    H = {}
    irofile = iter(a)
    for line in irofile:
        line = line.rstrip()
        if line.startswith('>'):
            prophage = genome+'_'+str(line[-1:])
            if prophage not in H:
                H[prophage] = {}
            SEQ = next(irofile)
            SEQ = SEQ[:-1]
            length = len(SEQ)
            total = length/3
            if length%3 == 0:
                counter = 1
                start = 0
                end = 3
                while counter <= total:
                    codon = SEQ[start:end]
                    start = start+3
                    end = end+3
                    if codon not in H[prophage]:
                        H[prophage][codon] = 1
                    elif codon in H[prophage]:
                        H[prophage][codon]+=1
                    counter +=1
            elif length%3 != 0:
                print 'error, not in multiples of 3'
    return H

def genome_profile(file_):
    a = open('/home3/hkang408/genome_ORFs/'+str(file_))
    genome = file_[:file_.find('.fasta')]
    H = {}
    irofile = iter(a)
    for line in irofile:
        line = line.rstrip()
        if line.startswith('>'):
            if genome not in H:
                H[genome] = {}
            SEQ = next(irofile)
            SEQ = SEQ[:-1]
            length = len(SEQ)
            total = length/3
            if length%3 == 0:
                counter = 1
                start = 0
                end = 3
                while counter <= total:
                    codon = SEQ[start:end]
                    start = start+3
                    end = end+3
                    if codon not in H[genome]:
                        H[genome][codon] = 1
                    elif codon in H[genome]:
                        H[genome][codon]+=1
                    counter +=1
            elif length%3 != 0:
                print 'error, not in multiples of 3'
    return H

a = open('/home3/hkang408/pp_fastas/tRNA_in_pp.txt', 'rU')

for line in a:
    i = line.split('\t')
    pp = i[0].rstrip()
    genome = pp[:pp.find('_')]
    fasta = genome+'.fasta'
    H_pp = codon_profile(fasta)   # for pp
    H_genome = genome_profile(fasta)    # for genome
    codon_pp = H_pp[pp][i[5].rstrip()]
    # SET NEW VALUES FOR HOST
    
    for prophage in H_pp:
        if prophage == pp:
            for codon in H_pp[prophage]:
                old_value = H_genome[genome][codon]
                new_value = int(old_value) - int(H_pp[prophage][codon])
                H_genome[genome][codon] = new_value

    codon_count_in_pp = H_pp[pp][i[5].rstrip()]
    avg_codon_in_pp = float(int(codon_count_in_pp)*1.0/sum(H_pp[pp].values()))
    codon_count_in_host = H_genome[genome][i[5].rstrip()]
    adjusted_codon_count_in_host = codon_count_in_host - 
    avg_codon_in_host = float(int(codon_count_in_host)*1.0/sum(H_genome[genome].values()))

    
    print str(codon_usage_in_pp) +': codon usage in pp'
    print str(codon_avg_in_pp) + ': avg codon usage in pp'
    print str(codon_usage_in_host) + ': codon usage in host'
    print str(codon_avg_in_host) + ': avg codon usage in host'
