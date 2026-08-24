#!/usr/bin/env python
import argparse

###############################################################################################################################################################

parser = argparse.ArgumentParser(description = 'Extract DHS associated with genes from a defined list')
parser.add_argument('genelist', type = str, help = 'Gene list in one column')
parser.add_argument('coords', type = str, help = 'Bed file with gene name in 4th column')
parser.add_argument('outfile', type = str, help = 'Output file')

args = parser.parse_args()

###############################################################################################################################################################

genelist = open(args.genelist, 'r').read()
genelist = genelist.splitlines()

outfile = open(args.outfile, 'w')
 
with open(args.coords, 'r') as coords:
    for x in coords:
        columns = x.strip().split('\t')       

        try:
            if columns[3] in genelist:
                outfile.write('\t'.join(columns[0:3]) + '\n')
        except IndexError:
            # Catch the error when there are fewer than 4 columns
            continue

outfile.close()

###############################################################################################################################################################