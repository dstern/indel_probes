from string import *

def get_in_oligo(sequence):
    seq = sequence.replace('-','')
    oligo = seq[20:80]
    GC = oligo.count('G') + oligo.count('C') 
    end = len(seq) - 80
#    print GC
    x = 21
    while x < end:
        test_oligo = seq[x:x+60]
        test_GC = test_oligo.count('G') + test_oligo.count('C')
#        print test_GC
        
        if GC < test_GC:
            GC = test_GC
            oligo = test_oligo
        
        x += 1
#    print GC
    if GC>30: #want between 23 and 30 GC
        return "skip"
    elif GC<22:
        return "skip"
    else:
        return oligo    

def get_del_oligo(sequence):
    seq = sequence.replace('-','') 
    oligo = seq[10:70]
    GC = oligo.count('G') + oligo.count('C') #want between 23 and 30 GC
    x = 11
    while x < 40:
        test_oligo = seq[x:x+60]
        test_GC = test_oligo.count('G') + test_oligo.count('C')
        
        if GC < test_GC:
            GC = test_GC
            oligo = test_oligo
        x += 1

    if GC>30:
        return "skip"
    elif GC<22:
        return "skip"
    else:
        return oligo    

def get_in_oligo_name(chromosome,start,species):
        
    if start < 100:
        position = '000000' + str(start)
    elif start < 1000:
        position = '00000' + str(start)
    elif start < 10000:
        position = '0000' + str(start)
    elif start < 100000:
        position = '000' + str(start)
    elif start < 1000000:
        position = '00' + str(start)
    elif start < 10000000:
        position = '0' + str(start)
    else:
        position = str(start)
        
    name = str(chromosome) + '_' + position + '_'+ species+ '_' + 'i'
    return name


def get_del_oligo_name(chromosome,start,species):
    if start < 100:
        position = '000000' + str(start)
    elif start < 1000:
        position = '00000' + str(start)
    elif start < 10000:
        position = '0000' + str(start)
    elif start < 100000:
        position = '000' + str(start)
    elif start < 1000000:
        position = '00' + str(start)
    elif start < 10000000:
        position = '0' + str(start)
    else:
        position = str(start)

    name = str(chromosome) + '_' + position + '_'+ species+ '_' + 'd'
    return name
