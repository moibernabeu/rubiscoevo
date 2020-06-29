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
parser.add_option("-a", "--analysys", dest="analysis",
                  help="Execute a model selection analysis and \
                  fittest model tree inference (y/ n).",
                  metavar="ANALYSIS")
(options, args) = parser.parse_args()

if options.file is None:
    file = input('Fasta or Phylip file containing sequences: ')
else:
    file = options.file

# Alignment
cmd_al = str("t_coffee -in " + file + " -method muscle_msa,mafft_msa,\
             clustalo_msa,clustalw2_msa -mode mcoffee")

print('Executed:' + cmd_al)
os.system(cmd_al)

# Convert Fasta to Phylip if input is Fasta
file_name = str(file.split('.')[0])

aln_name = str(file_name + '.aln')
aln_ext = str(aln_name.split('.')[-1])

if aln_ext == 'aln':
    cmd_readal = str('readal -in ' + aln_name + ' -out ' +
                     file_name + '_aligned.fa -fasta')
    os.system(cmd_readal)
pass

# aln_name_phy = str(file_name + '_aligned.phy')
if options.seqtype is None:
    type = input('Data type of your alignment (DNA/ AA): ')
else:
    type = options.seqtype

dir = os.getcwd()

cmd_trim = str('java -jar \
               /Users/luca/software/BMGE-1.12/BMGE.jar -t ' +
               type + ' -i ' + dir + '/' + file_name +
               '_aligned.fa -op ' + dir + '/' +
               file_name + '_trimmed.phy')
print('Executed: ' + cmd_trim)
os.system(cmd_trim)

os.system('rm *.dnd')
os.system('rm *.html')
os.system('rm *.aln')

# Continue with model selection implemented in IQ-TREE
if options.analysis is None:
    sel = input('Do you want to continue with model\
    selection analysis? (y/ n): ')
else:
    sel = options.analysis

if sel == "y":
    if type == 'DNA':
        cwd = os.getcwd()

        file = str(file_name + '_trimmed.phy')

        if options.outgroup is None:
            outgroup = input('Define outgroup sequence: ')
        else:
            outgroup = options.outgroup

        iqtree = str('iqtree -s ' + file + ' -bb 2000 -st DNA -m TEST \
                     -t BIONJ -o ' +
                     outgroup + ' -mem 4G -wt -pre model_sel')
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

        file = str(file_name + '_trimmed.phy')

        if options.outgroup is None:
            outgroup = input('Define outgroup sequence: ')
        else:
            outgroup = options.outgroup

        iqtree = str('iqtree -s ' + file + ' -mset Blosum62,cpREV,\
                     Dayhoff,DCMut,JTT,JTTDCMut,Poisson,PMB,VT,WAG,\
                     GTR20 -t BIONJ -o ' +
                     outgroup + ' -mem 4G -wt -pre model_sel')
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
