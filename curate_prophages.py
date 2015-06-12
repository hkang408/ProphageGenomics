def pp_get(path,name):

# Setting all necessary parameters
    count = 0
    pp_counter = 0
    switch = 0 # pp designator
    int_pos = 0
    hypo_counter = 0
    t_id = ''
    t_name = ''         # t for temp, l for left, r for right

    l_id = ''
    l_name = ''
    b = open('pp_summary.txt', 'a+')
    try:
        a = open(path+'/prophage_tbl.txt', 'rU')
    except:
        print str(path) + ': file not found'
        return
    
    irofile = iter(a)
    int_id = 0
    int_pos = 0
    repeat_l = 0
    repeat_r = 0
    repeat_reason = 0
    
    for line in irofile:
        i = line.split('\t')
        
        if switch == 0: # not in pp
            if i[9].startswith('0'): ## Is not a phage
                pass
            elif i[9].startswith('1'): ## Beginning of phage
                switch = 1
                pp_counter += 1
                if 'hypothetical' in i[1] or 'putative' in i[1] or 'Hypothetical' in i[1] or 'Putative' in i[1]:
                    hypo_counter += 1
                if int(i[3]) < int(i[4]):
                    start = i[3]
                elif int(i[3]) > int(i[4]):
                    start = i[4]
                else:
                     print 'error in setting start position'
                         
                contig = i[2]
                try:                    # This only works when the repeat designation is at the start of pp region...
                    repeat_l = i[14]
                    repeat_r = i[15]
                    repeat_reason = i[16]
                    repeat_reason = repeat_reason.rstrip()
                except:
                    repeat_l = 0
                    repeat_r = 0
                    repeat_reason = 0
                if ('integrase' in i[1] or 'Integrase' in i[1]):        # Not all will have an integrase... and some may have multiple.
                    int_pos = pp_counter
                    int_id = i[0]

                    
        elif switch == 1:
            if (i[9].startswith('1') and contig == i[2]): #Same prophage continued.
                pp_counter +=1
                if 'hypothetical' in i[1] or 'putative' in i[1] or 'Hypothetical' in i[1] or 'Putative' in i[1]:
                    hypo_counter += 1
                if int(i[3]) > int(i[4]):
                    end = i[3]
                elif int(i[3]) < int(i[4]):
                    end = i[4]
                else:
                     print 'error in setting end position'

                
                if ('integrase' in i[1] or 'Integrase' in i[1]):
                    int_pos = pp_counter                # This would overwrite any pp found at beginning of integrase
                    int_id = i[0]               

            elif (i[9].startswith('1') and contig != i[2]):   # contigs changed in pp. .: 2 prophages
                count +=1
                try:
                    r_id = 0
                    r_name = 0
                    pp_len = int(end) - int(start) # Uses previous end value
                    total_genes = pp_counter
                    integrase_position = (int(int_pos)*1.0)/int(pp_counter)*1.0
                    hypo_perc = (int(hypo_counter)*1.0)/int(total_genes)*1.0
                except:
                    pass
                
                if int(pp_counter) < 5:
                    print str(i[0]) + ': Region too small.'
                    int_id = 0
                    int_pos = 0
                    repeat_l = 0
                    repeat_r = 0
                    repeat_reason = 0
                    pp_counter = 0
                    start = 0
                    hypo_counter = 0
                    count -= 1
                elif int(pp_counter) >= 5:
                    b.write(
                    name+'*' + '\t' + str(count) + '\t' + str(contig) + '\t' + str(start) + '\t' + str(end) + '\t'
                    + str(l_id) + '\t' + str(l_name) + '\t' + str(r_id) + '\t' +str(r_name) + '\t' + str(hypo_perc) + 
                    '\t')

                    if int_id != 0:                                     # If integrase has been found
                        b.write(str(int_id) + '\t' + str(integrase_position) + '\t')
                    else:
                        b.write('\t\t')

                    if repeat_l !=0:                                    # If successful in finding repeat region
                        b.write(str(repeat_l) + '\t' + str(repeat_r) + '\t' + str(repeat_reason))
                    else:
                        b.write('\t\t')
                    b.write('\n')
                
                int_id = 0
                int_pos = 0
                repeat_l = 0
                repeat_r = 0
                repeat_reason = 0
                pp_counter = 1
                start = 0
                hypo_counter = 0
                contig = i[2]
                if int(i[3]) < int(i[4]):
                    start = i[3]
                elif int(i[3]) > int(i[4]):
                    start = i[4]
                else:
                    print 'error in setting start position'
                
            elif i[9].startswith('0'): # Reached the end of pp
                count+=1
                try:
                    r_id = i[0]
                    r_name = i[1]
                    pp_len = int(end) - int(start)
                    total_genes = pp_counter
                    integrase_position = (int(int_pos)*1.0/int(pp_counter))*1.0
                    hypo_perc = (int(hypo_counter)*1.0)/int(total_genes)*1.0
                except:
                    pass
                
                if int(pp_counter) < 5:
                    int_id = 0
                    int_pos = 0
                    repeat_l = 0
                    repeat_r = 0
                    repeat_reason = 0
                    pp_counter = 0
                    start = 0
                    hypo_counter = 0
                    count -= 1
                    switch = 0

                elif int(pp_counter) >= 5:
                    skew = get_AT(name,contig,start,end)                    
                        
                    b.write(
                    name + '\t' + str(count) + '\t' +str(contig) + '\t' + str(pp_len) + '\t'
                    + str(hypo_perc)+'\t'+ str(skew) + '\t')
                    
                    if int_id != 0:                                     # If integrase has been found
                        b.write(str(int_id) + '\t' + str(integrase_position) + '\t')
                    else:
                        b.write('\t\t')

                    if repeat_l !=0:                                    # If successful in finding repeat region
                        b.write(str(repeat_l) + '\t' + str(repeat_r) + '\t' + str(repeat_reason))
                    else:
                        b.write('\t\t')

                    b.write('\n')
                
                int_id = 0
                int_pos = 0
                repeat_l = 0
                repeat_r = 0
                repeat_reason = 0
                pp_counter = 0
                start = 0
                hypo_counter = 0
                switch = 0 #No longer in prophage
    a.close()
    b.close()
    return


def automate():
    a = open('complete_genomes.txt','U')
    b = open('pp_summary.txt', 'w')
    b.write('Name\tpp count\tContig\tLength\tHypo Perc\tAT Perc\tIntegase ID\tIntegrase Position\tLeft Repeat\tRight Repeat\tRepeat Reason\n') 
    b.close()
    for i in a:
        i = i.strip()
        i = i.replace('\n','')
        path = '/home3/katelyn/KBase/prophages/SEED/'+i
        pp_get(path,i)
        
    a.close()
    return


def get_AT(genome,contig,start,end):
    a = open('/home3/katelyn/KBase/prophages/SEED/'+str(genome)+'/contigs')
    irofile = iter(a)
    for line in a:
        if '>'+str(contig) == line.rstrip():
            sequence = next(irofile)
            SEQ = sequence[int(start):int(end)]
            T = SEQ.count('T')
            A = SEQ.count('A')
            skew = (float(A)+float(T))/len(SEQ)

    return skew

automate()
