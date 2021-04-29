"""
Run this with python folding_with_reactivities.py.  This is written for Python2, but also works in Python3
This program runs the folding of RNA sequences with SHAPE data.  SHAPE reactivities are incorporated as an energetic contribution to the folding energy following Deigan et al. [PNAS 106 97 (2009)]
The intercept (b) and slope (m) parameterising the contribution from the shape reactivities are read from the list given in m_and_b.dat 
The shape reactivities are supplied in a single column file (reactivities_file.shp), with the reactivity on the nth line representing the reactivity of the nth nucleotide in the sequence.
A reactivity should be provided for every nucleotide, where nucleotides without available reactivties are indicated with a -999.
The sequence should be provided in a fasta file, with the sequence name read from "> | | | sequence_name"
An output directory "extracted" should be created before running
RJB
"""

import subprocess
import sys

m_and_b_file = 'm_and_b.dat' 

genome_file = 'sequence.fasta'

fold_sample_size = 1000         

shp_file = 'reactivities_file.shp'

m_and_b_list = [[j for j in i.strip('\n').split(' ') if j!=''] for i in open(m_and_b_file).readlines()]

m_and_b_list = [j[0]+' '+j[1] for j in m_and_b_list]

seq_name = open(genome_file).readline()[1:].split('|')[3].replace(' ','').strip('\n')

for m_and_b in m_and_b_list:
    m = m_and_b.split(' ')[0]
    b = m_and_b.split(' ')[1]
    print('folding '+m_and_b)
    file_name = 'extracted/'+seq_name+'.folds'+str(m)+'_'+str(b)

    bashCommand = './Sfold_mod.x -p %d -i %s -o %s -s %d -shp %s -m %s -b %s ' % (fold_sample_size, genome_file, file_name, 123456789,shp_file, m ,b)

    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    process.communicate()

