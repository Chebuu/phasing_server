#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2019-Sep-25, latest modified on 2019-Sep-25
# Apply a uniform bfactor to all residues

# Example in Linux: python apply_uniform_bfactor.py a.pdb 1.5 1_uniform.pdb

import numpy as np
import sys
import argparse


def apply_bfactor(topology, bfactor, output):

    if topology[-4:].lower() != ".pdb":
        topology = topology + ".pdb"
    topology_prefix = topology.split(".")[0]
    bfactor = float(bfactor)
    # bfactor_text = str('{:.2f}'.format(float(bfactor)))
    # print(bfactor_text)
    with open ('%s' %topology, 'r') as fopen:
        lines = fopen.readlines()
        data = []
        for line in lines:
            if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
                if line[12:16].strip() != "HB":                
                    data.append("%s%6.2F%s" % (line[:60], bfactor, line[66:]))

    with open('%s' %output, 'w') as f2:
        for line in data:
            f2.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="This script applys an uniform B-factor value to the input structure. \
        An example for this script is python apply_uniform_bfactor.py best.pdb 1.5 best_uniform.pdb")
    parser.add_argument("topology", help="The file name of input structure", type=str)
    parser.add_argument("bfactor", help="The value of B-factor", type=float)
    parser.add_argument("output", help="The output file name", type=str)


    args = parser.parse_args()
    topology = args.topology
    bfactor = args.bfactor
    output = args.output
    apply_bfactor(topology, bfactor, output)

if __name__ == '__main__':
    main()
