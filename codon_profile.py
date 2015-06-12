### Input takes in the modal codon usage program http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2839124/
### Identifies the codon usage of codon and amino acid in question, returns modal usage freq.

def make_tmp_pp(prophage):
    prophage = prophage.rstrip()
    genome = prophage[:prophage.find('_')]
    pp = prophage[prophage.find('_')+1:]
#    a = open('/home3/hkang408/pp_fastas/'+ genome+'.fasta', 'rU')
    a = open('/home3/hkang408/pp_genes/'+genome+'.fasta', 'rU')
    b = open('/home3/hkang408/pp_codon_profile/'+prophage+'.fasta', 'w')
    irofile = iter(a)
    for line in irofile:
        line = line.rstrip()
        if line.startswith('>'):
            prophage = line[line.find('_')+1:]
            if prophage == pp:
                b.write(line+'\n'+next(a))
    a.close()
    b.close()

def codon_output_analyzer(file_,AA_,CODON):
    codons = ['GCA', 'GCG', 'GCT', 'GCC', 'TGT', 'TGC', 'GAT', 'GAC', 'GAA', 'GAG',
              'TTT', 'TTC', 'GGA', 'GGG', 'GGT', 'GGC', 'CAT', 'CAC', 'ATA', 'ATT', 'ATC', 'AAA', 'AAG',
              'TTA', 'TTG', 'CTA', 'CTG', 'CTT', 'CTC', 'AAT', 'AAC', 'CCA', 'CCG', 'CCT', 'CCC', 'CAA', 'CAG', 'AGA', 'AGG', 'CGA', 'CGG', 'CGT', 'CGC',
              'AGT', 'AGC', 'TCA', 'TCG', 'TCT', 'TCC', 'ACA', 'ACG', 'ACT', 'ACC', 'GTA', 'GTG', 'GTT', 'GTC', 'TAT', 'TAC', 'ATG', 'TGG']
    AA = ['Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile', 'Lys', 'Leu', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr', 'Val', 'Tyr', 'Met', 'Trp']

    a = open(file_,'rU')
    out = a.readline()
    out = out.split('\t')
    i = out[1].split('|')
    counter1 = 0
    counter2 = 0

    H = {}
    for aa in i:
        aa = aa.rstrip()
        H[AA[counter1]] = {}
        codon = aa.split(',')
        for each in codon:
            H[AA[counter1]][codons[counter2]] = each
            counter2 +=1
        counter1 +=1
    return H[AA_][CODON]

def tRNA():
    a = open('tRNA_in_pp.txt', 'rU')
    b = open('codon_differences.txt', 'w')
    b.write('Prophage\tCodon\tPhage Usage\tHost Usage\n')
    convert_codon = {'G':'C', 'A':'T', 'T':'A', 'C':'G'}
    for line in a:
        i = line.split('\t')
        pp = i[0].rstrip()
        aa = i[4].rstrip()
        codon = i[5].rstrip()
        genome = pp[:pp.find('_')]
        try:
            op_codon = convert_codon[codon[2]]+convert_codon[codon[1]]+convert_codon[codon[0]]
            host_codon = codon_output_analyzer('/home3/hkang408/codon_usage_programs/'+genome+'.modal_freq', aa, op_codon)
            phage_codon = codon_output_analyzer('/home3/hkang408/pp_codon_profile/'+pp+'.modal_freq', aa, op_codon)
            if float(phage_codon) > float(host_codon):
                call = '+'
            else:
                call = '-'

            b.write(pp+'\t'+op_codon+'\t'+str(phage_codon)+'\t'+str(host_codon)+'\t'+str(call)+'\n')
            
        except:
            #print aa +':\tnot an accepted AA'
            pass
        





'''

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

'''
