"""
This script is run as python stem_loop_motif_frequency.py
It is written for python2, but should work in python3, provided the required modules are installed.
This script reads in the sampled fold files before searching for stemloops that overlap with the sequence motifs that are supplied in the 'motif_list' file.
The motifs for searching are supplied in the 'motif_list' file, as pos sequence pairs.
The output is the frequency of the motif occuring in stemloop loop sections (out of a maximum of 1000 - based on the number of samples requested in the earlier script).
This is written in the MT_pos_sequence file.  The picked_sls_pos_sequence file records the particular stemloops included. 
RJB

"""

import numpy as np
import re
import sys
from subprocess import PIPE, Popen

def database_maker(filename,seq):
    
    db = [{} for pos in range(len(seq))]
    
    folds = [i.strip('\n').split(' ')[0] for i in open(filename,'r').readlines()]

    min_structure = '\(\.+?\)'
    
    for fold in folds:
        for s_m in re.finditer(min_structure,fold):
            motif = seq[s_m.start()+1:s_m.end()-1]
            pos = s_m.start()+1
            if motif in db[pos].keys():
                db[pos][motif] += 1
            else:
                db[pos][motif] = 1

    return db

def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0]

if __name__ =='__main__':
    
    file_list = [ i for i in cmdline('ls ./extracted/*folds*').split('\n') if i!='' and 'output' not in i ]

    seq = open(file_list[0],'r').readline().strip('\n')

    m_and_b_list =  [[ float(f.split('_')[-2].split('folds')[-1]) , float(f.split('_')[-1]) ] for f in file_list ]

    dmb = [[ float(f.split('_')[-2].split('folds')[-1]) , float(f.split('_')[-1]) , database_maker(f,seq)] for f in file_list ] 

    all_folds = [ {} for i in range(len(seq)) ] 
    

    for db in dmb:
        for pos in range(len(db[-1])):
            for seq in db[-1][pos].keys():
                all_folds[pos][seq] = all_folds[pos].get(seq,[]) + [[db[0],db[1],db[-1][pos][seq]]]
                    
    tmp_list = []
    for i in range(len(dmb)):
        t1=[]
        for j in range(len(dmb[i][-1])):
            t1.extend([dmb[i][-1][j][seq] for seq in dmb[i][-1][j].keys()])
        tmp_list.append([dmb[i][0],dmb[i][1],np.mean(t1)])
    
    m_and_b_list = tmp_list
    
    f=open('motif_list','r')
    motifs = [[i.split('\t')[-1].strip('\n'),int(i.split('\t')[0])-1] for i in f.readlines()]
    f.close()

    d_val={}
    for val in m_and_b_list:  
        d_val[str(val[0])+'_'+str(val[1])]=[0 for m in motifs]
    for motif in motifs:
        d1=0.0
        f=open('MT_'+str(motif[1]+1)+'_'+motif[0],'w')
        f.write('m,b,freq\n')
        motif_pos = set([motif[1]+i for i in range(len(motif[0]))])
        picked_sls = []
        for pos in range(motif[1]-25,motif[1]+len(motif[0])+1):
           for seq in all_folds[pos].keys():
               seq_pos = set([pos+i for i in range(len(seq))])
               if len(seq_pos.intersection(motif_pos))!=0:
                   picked_sls.append([pos,seq,max([i for i in all_folds[pos][seq]],key=lambda x: x[-1]),len(all_folds[pos][seq]),np.mean([i[-1] for i in all_folds[pos][seq]])])
        g=open('picked_sls_'+str(motif[1]+1)+'_'+motif[0],'w')
        for pick in picked_sls:
            g.write(str(pick[0])+'\t'+str(pick[1])+'\t'+str(round(pick[4],2))+'\t'+str(pick[3])+'\t'+str(pick[2])+'\n')
        g.close()
        motif_d ={}
        for sl in picked_sls:
            for data in all_folds[sl[0]][sl[1]]:
                if str(data[0])+'_'+str(data[1]) in motif_d.keys():
                    motif_d[str(data[0])+'_'+str(data[1])] += data[2]
                else:
                    motif_d[str(data[0])+'_'+str(data[1])] = data[2]
                    
        for val in m_and_b_list:
            written=False
            if str(val[0])+'_'+str(val[1]) in motif_d.keys():
                data = [val[0],val[1],motif_d[str(val[0])+'_'+str(val[1])]]
                if d1!=data[0]:
                   #f.write('\n')   # Uncommenting these lines will add line breaks between m values, to aid in plotting.
                   f.write(str(data[0])+','+str(data[1])+','+str(data[2])+'\n')
                   d1 = data[0]
                   written=True 
                else:
                   f.write(str(data[0])+','+str(data[1])+','+str(data[2])+'\n')
                   written=True 
            if not written:
                if d1!=val[0]:
                   #f.write('\n')
                   f.write(str(val[0])+','+str(val[1])+','+str(0)+'\n')
                   d1 = val[0]
                   written=True 
                else:
                   f.write(str(val[0])+','+str(val[1])+','+str(0)+'\n')
                   written=True 
        f.close()
        
