#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013 Junko Tsuji

# IO functions of heatlogoSS.

import sys, commands
from StringIO import StringIO

from weblogolib import read_seq_data
from corebio.utils import find_command
from corebio.seq import unambiguous_protein_alphabet


# read protein sequence and make Seq class
def readSequence(f, ip):
    if f == sys.stdin:
        fin = StringIO(sys.stdin.read())
    else:
       fin = open(f, "r")
    seqs = read_seq_data(fin, ip)
    if seqs.alphabet != unambiguous_protein_alphabet:
        raise Exception("input sequences should be protein seqeunces")
    return seqs


# generate single fasta files
def makeSingleFasta(seqs, wd):
    L = len(seqs)
    linesize = 50
    files, exclude = [], []
    for i in range(L):
        name = "seq" + str(i)
        seq = [s for s in seqs[i].upper() if s in "ACDEFGHIKLMNPQRSTVWY"]
        seq = "".join(seq)
        seqlen = len(seq)

        if seqlen <= 20:
            exclude.append(i)
            continue

        file = wd + "/" + name + ".fa"
        f = open(file, "w")
        f.write(">" + name + "\n")
        beg = 0
        while beg < seqlen:
            end = beg + linesize
            f.write(seq[beg:end])
            beg = end
        f.close()

        files.append(file)
    return files, exclude

