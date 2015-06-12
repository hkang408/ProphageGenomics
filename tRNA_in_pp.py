#   This script is used to identify tRNA sequences in the prophage. 
#



def _query(argA):
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    pipe = Popen(['/home3/katelyn/opt/tRNAscan-SE/tRNAscan-SE', '-P',  argA
   ],             #        '-dust', 'no'
    stdin=PIPE,
    stdout=PIPE)
    print pipe.communicate()[0]

def identify_tRNAs():            # python identify_tRNAs >out_tRNA
    a = open('pp_fastas.txt','rU')

    path = '/home3/hkang408/pp_fastas/'

    for line in a:
        line = line.rstrip()
        argA = path+str(line)
        _query(argA)

    a.close()

def tRNA_data_sum():       #   Run after: grep 'g\.' out_tRNA >tRNA_in_pp.txt
    a = open('tRNA_in_pp.txt', 'rU')
    H = {}

    for line in a:
        line = line.rstrip()
        i = line.split('\t')
        pp = i[0]
        start = i[2]
        end = i[3]
        AA = i[4]
        score = float(i[8])
        if pp not in H:
            H[pp] = [AA]
        elif pp in H:
            H[pp].append(AA)

    print len(H)    # Number of unique phages w/ tRNA

    list_ = []
    for each in H:
        for j in H[each]:
            list_.append(j)

    X = {}
    for i in list_:
        if i not in X:
            X[i] = str(list_.count(i))

    for each in X:                  # AA Distribution
        print each + '\t' + X[each]
