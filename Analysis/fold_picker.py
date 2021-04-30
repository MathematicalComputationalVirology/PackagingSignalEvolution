"""
This script runs using python fold_picker.py folds_file
This script is written to run in python2 and requires the selected folds_file (output from an earlier script) to be given on the command line
This script returns the minimum energy fold associated with a given m & b combination.  The fold is output in dot-bracket format  
It also returns the number of base pairs that this structure shares with a consensus_fold generated from the sampled folds 
RJB
"""

#python fold_picker.py folds_file 
import sys
import os

def pp_list(fold):
    pp2 = []
    pp = [-1 for i in range(len(fold))]
    stack = [0 for i in range(len(fold))]
    t1 = 0
    for i in range(len(fold)):
        if fold[i]=='(':
            t1 += 1
            stack[t1] = i
        elif fold[i] ==')':
            j = stack[t1]
            pp[i] = j
            pp[j] = i
            t1 -= 1 
    for pos in range(len(fold)):
        if pp[pos]>=0:
            pp2.append([pos+1,pp[pos]+1])
        else:
            pp2.append([pos+1,pp[pos]])
    pp2 = [x if x[-1]!=-1 else [x[0],0] for x in pp2]
    return pp2 

f=open(sys.argv[1],'r')
seq = f.readline().strip('\n')
folds_and_energies = [[i.split(' ')[0],float(i.split(' ')[-1])] for i in f.readlines()]
folds = [i[0] for i in folds_and_energies]
energies = [i[1] for i in folds_and_energies]
f.close()

m=float(sys.argv[1].split('folds')[-1].split('_')[0])
b=float(sys.argv[1].split('folds')[-1].split('_')[1])

bp_dict={}

for i in range(len(folds_and_energies)):
    bps=pp_list(folds_and_energies[i][0])
    bps = [sorted(j) for j in bps]
    bp_keys = list(set([str(j[0])+'_'+str(j[1]) for j in bps]))
    for bp in bp_keys:
            if bp in bp_dict.keys():
                    bp_dict[bp] += 1
            else:
                    bp_dict[bp] = 1 
    folds_and_energies[i].append(pp_list(folds_and_energies[i][0]))

fold=['p' for i in range(len(folds[0]))]

for bp in bp_dict.keys():
        bp_1 = int(bp.split('_')[0])
        bp_2 = int(bp.split('_')[1])
        if bp_dict[bp] >= int(0.5*len(folds)):
                if bp_1==0:
                        if fold[bp_2-1]!='p':
                                print 'error',bp_1,bp_2
                        fold[bp_2-1]='.'
                else:
                        if fold[bp_2-1]!='p':
                                print 'error2',bp_1,bp_2
                        fold[bp_1-1]='('
                        fold[bp_2-1]=')'

consensus_fold_0 =''.join(fold)
consensus_fold = consensus_fold_0.replace('p','.')

con_bp = pp_list(consensus_fold)
con_bp = [sorted(j) for j in con_bp]
con_bp_set = set([str(i[0])+'_'+str(i[1]) for i in con_bp if i[0]!=0])

for l in folds_and_energies:
        l1 = [sorted(j) for j in l[2]]
        l1 = set([str(j[0])+'_'+str(j[1]) for j in l[2] if j[1]!=0])
        l.append(len(l1.intersection(con_bp_set)))

folds_and_energies.sort(key=lambda x:x[1])

min_sampled_fold = folds_and_energies[0][0]
min_sampled_energy = folds_and_energies[0][1]
min_sampled_bp_con = folds_and_energies[0][-1]

f=open(sys.argv[1]+'_output_folds','w')
f.write('>min_sampled_fold '+str(min_sampled_energy)+' '+str(min_sampled_bp_con)+'\n')
f.write(seq+'\n')
f.write(min_sampled_fold+'\n')
