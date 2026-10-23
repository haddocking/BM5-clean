#!/usr/bin/env python

"""
Simplify izone files to speed up profit usage

# Usage: python join_zones.py <'XYZ'_refe.izone>
# Max. number of residues per chain: 9999

Author: {0} ({1})
Initial work: Leander Goldbach (l.d.goldbach@students.uu.nl)
Modified (feb. 2017): Mikael Trellet (mikael.trellet@gmail.com)

This program is part of the analysis tools distributed with HADDOCK
or with the HADDOCK tutorial.
"""

import sys
import os
import re

__author__ = "Leander Goldbach, Mikael Trellet"
__email__ = "l.d.goldbach@students.uu.nl, mikael.trellet@gmail.com"

USAGE = __doc__.format(__author__, __email__)


def makelist(inputfile):
    """ take .izone file and creates sorted list of all residues for chain A and B """
    lists = {}
    prev_chain = None
    chains = []
    for line in inputfile:
        chain = line.split()[1][0]
        if chain != prev_chain:
            chains.append(chain)
            lists[chain] = []
        resid = re.findall(r"\d+", line)[0]
        lists[chain].append(int(resid))
        prev_chain = chain
    for ch in chains:
        lists[ch].append(ch)
    for l in lists:
        # Python 2 ordering: residue numbers first, then the chain letter
        lists[l] = sorted(x for x in lists[l] if isinstance(x, int)) + [x for x in lists[l] if isinstance(x, str)]
    # lists[chains[0]].sort()
    # lists[chains[1]].sort()
    return lists



def combine_zones(listX):
    """takes sorted lists and combines neighbouring residues,
    doubles individual residues, e.g. [2, 3, 4, 12] -> [2, 4, 12, 12]"""
    Xzones=[listX[0]]
    iden = listX[len(listX)-1]
    del listX[len(listX)-1]
    x = 1
    for num in listX:
        if  x == len(listX):
            break
        elif num + 1 == listX[x]:
            x = x + 1
        else:
            Xzones.append(listX[x-1])
            Xzones.append(listX[x])
            x = x + 1
    Xzones.append(listX[len(listX) - 1])
    return Xzones



def output(zones):
    """creates a file from both lists, two list objects become one zone, e.g.
    [2, 4, 12, 12] -> 'zone A2-A4 \n zone A12-A12 \n'"""
    c = 65
    for z in zones:
        for position, value in enumerate(z):
            if position == 0 or position%2 == 0:
                line = "zone {}.{}-".format(chr(c), str(value))
            elif position%2 != 0:
                line += "{}.{}\n".format(chr(c), str(value))
                yield line
        c += 1
    # for position, value in enumerate(zonesB):
    #     #print position, value
    #     if position == 0 or position%2 == 0:
    #         line = "zone B.{}-".format(str(value))
    #     elif position%2 != 0:
    #         line += "B.{}\n".format(str(value))
    #         yield line


def check_input(args):
    """Checks whether to read from stdin/file and validates user input/options."""

    if len(args) < 1:
        # No file, from pipe
        if not sys.stdin.isatty():
            izonefh = sys.stdin
        else:
            sys.stderr.write(USAGE)
            sys.exit(1)
    elif len(args) == 1:
        # File
        if not os.path.isfile(args[0]):
            sys.stderr.write('File not found: ' + args[0] + '\n')
            sys.stderr.write(USAGE)
            sys.exit(1)
        else:
            izonefh=open(args[0], 'r')
    else:
        sys.stderr.write(USAGE)
        sys.exit(1)

    return izonefh

if __name__ == '__main__':
    # Check Input (file or pipe)
    izonefh = check_input(sys.argv[1:])
    
    # Create sorted residue lists
    res_lists = makelist(izonefh)
    #print Alist, Blist

    # Create combined zones for continuous ranges of resids
    zones = []
    for l in res_lists:
        zones.append(combine_zones(res_lists[l]))
    # Azones = combine_zones(Alist)
    # Bzones = combine_zones(Blist)

    # Generate output
    new_izone = output(zones)

    try:
        sys.stdout.write(''.join(new_izone))
        sys.stdout.flush()
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

    izonefh.close()

