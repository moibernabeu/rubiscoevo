#!/usr/bin/env python3

# Moisès Bernabeu
# Started Burjassot February, 2020
# Phylogenetic analysis pipeline
# Reuieres: t_coffee (muscle, mafft, clustalo, clustalw2), readal, BMGE
# and IQ-TREEp

import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="file",
                  help="Input fasta (.fasta, .fa) file", metavar="FILE")
parser.add_option("-t", "--seq_type", dest="seqtype",
                  help="Name of sequence to filter (DNA/ AA)",
                  metavar="SEQUENCE TYPE")
parser.add_option("-o", "--outgroup", dest="outgroup",
                  help="Outgroup to phylogenetic analysis",
                  metavar="OUTGROUP")
(options, args) = parser.parse_args()

if options.file is None:
    file = input('Fasta or Phylip file containing sequences: ')
else:
    file = options.file

# Continue with model selection implemented in IQ-TREE
type = options.seqtype
if type == 'DNA':
    cwd = os.getcwd()

    file = str(file)

    if options.outgroup is None:
        outgroup = input('Define outgroup sequence: ')
    else:
        outgroup = options.outgroup

    iqtree = str('iqtree -s ' + file + ' -bb 2000 -st DNA -m TEST \
                 -t BIONJ -o ' +
                 outgroup + ' -mtree -wt -pre model_sel')
    print(iqtree)
    os.system(iqtree)
    os.system('mkdir selection')
    os.system('mv model_sel.* selection')
    os.chdir('selection')
    os.system('gunzip *.gz')
    os.system('mkdir txt')
    os.system('cp model_sel.* txt')
    os.chdir('txt')
    os.system('rename -a .txt *')

    output = open('model_sel.log.txt', 'r')
    output_read = output.read()
    table = open('selection_table.csv', 'w')
    splitted = str(output_read.split(') ...\n ')[-1])
    splitted = str(splitted.split('chosen according to BIC\n')[0])
    csv = splitted.replace('    ', ';')
    csv = csv.replace('   ', ';')
    csv = csv.replace('\n  ', '\n')
    csv = csv.replace('\n ', '\n')
    csv = csv.replace('  ', ';')
    csv = csv.replace(' ', ';')
    csv = csv.replace(';;', ';')
    csv = csv.replace(';;', ';')
    csv = csv.replace(';', ',')
    print(csv)
    table.write(csv)
    output.close()
    table.close()
elif type == 'AA':
    cwd = os.getcwd()

    file = str(file)

    if options.outgroup is None:
        outgroup = input('Define outgroup sequence: ')
    else:
        outgroup = options.outgroup

    iqtree = str('iqtree -s ' + file + ' -mset Blosum62,cpREV,Dayhoff,DCMut,JTT,JTTDCMut,Poisson,PMB,VT,WAG -t BIONJ -o ' +
                 outgroup + ' -mtree -bb 2000 -wt -pre model_sel')
    print(iqtree)

    os.system(iqtree)
    os.system('mkdir selection')
    os.system('mv model_sel.* selection')
    os.chdir('selection')
    os.system('gunzip *.gz')
    os.system('mkdir txt')
    os.system('cp model_sel.* txt')
    os.chdir('txt')
    os.system('rename -a .txt *')

    output = open('model_sel.log.txt', 'r')
    output_read = output.read()
    table = open('selection_table.csv', 'w')
    splitted = str(output_read.split(') ...\n ')[-1])
    splitted = str(splitted.split('chosen according to BIC\n')[0])
    csv = splitted.replace('    ', ';')
    csv = csv.replace('   ', ';')
    csv = csv.replace('\n  ', '\n')
    csv = csv.replace('\n ', '\n')
    csv = csv.replace('  ', ';')
    csv = csv.replace(' ', ';')
    csv = csv.replace(';;', ';')
    csv = csv.replace(';;', ';')
    csv = csv.replace(';;', ';')
    csv = csv.replace(';;', ';')
    csv = csv.replace(';;', ';')
    csv = csv.replace(';', ',')
    print(csv)
    table.write(csv)
    output.close()
    table.close()
else:
    print('Finished!')
