#!/usr/bin/env python3

import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', '--file', dest='file',
                  help='Input fasta (.fasta, .fa) file', metavar='FILE')
parser.add_option('-p', '--prefix', dest='prefix',
                  help='Prefix to the folder', metavar='PREFIX')
(options, args) = parser.parse_args()

if os.path.exists(str(options.prefix + '_model_trees/')) is False:
    os.mkdir(str(options.prefix + '_model_trees/'))
    pass
os.chdir(str(options.prefix + '_model_trees/'))

with open(str('../' + options.file), 'r') as file:
    for line in file:
        headers = line.split(':')[0]
        prototree = line.split(';\n')[0]
        tree = prototree.split(': ', 2)[-1]
        nfile = open(str(headers + '_tree.nwk'), 'w')
        nfile.write(tree + ';')
        nfile.close()
        pass
    pass

os.system(str('cp ../' + options.file + ' ../' + options.prefix +
          '_model_trees/'))
os.system(str('mv ' + options.file + ' 00_' + options.file))
