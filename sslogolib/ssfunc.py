#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2013 Junko Tsuji

# Functions to convert disorder and secondary structure results
# to Sequence Logo like format.


import commands, string
from numpy import array

# insert gaps to prediction symbols
def gapify(prob, ref):
    pos = 0
    L = len(ref)
    gapped = []
    for i in range(L):
        if ref[i].upper() not in "ACDEFGHIKLMNPQRSTVWY":
            gapped.append(array([0.0, 0.0, 0.0, 0.0]))
        else:
            gapped.append(prob[pos])
            pos += 1
    return gapped


# execute predictors
def exePredict(init, exe, fasta, outDat, ignore, index, clean):
    commands.getoutput(init)
    prob = []
    L = len(fasta)
    for i in range(L):
        p = []
        commands.getoutput(exe % fasta[i])
        outProb = commands.getoutput(outDat)
        # pack prediction probabilities
        for line in open(outProb, "r"):
            data = line.split()
            if line.startswith(ignore):
                pass
            elif not data:
                pass
            else:
                d = [ float(data[j]) for j in index ]
                p.append(d)
        prob.append(p)
        commands.getoutput(clean)
    return prob


# merge Poodle-L and PsiPred results:
def mergeResults(seqs, poodleProb, psipredProb):
    A = len(seqs)
    L = len(seqs[0])
    stProb = [[0,]*4 for i in range(0, L)]
    for i in range(len(poodleProb)):
        line = []
        for j in range(len(poodleProb[i])):
            psipredProb[i][j] = [p * 1.5 for p in psipredProb[i][j]]
            colm = poodleProb[i][j] + psipredProb[i][j]
            total = sum(colm)
            line.append(array([value/total for value in colm]))
        gapped = gapify(line, seqs[i])
        for k in range(0, L):
            stProb[k] += gapped[k]
    for i in range(0, L):
        stProb[i] /= float(A)
    return stProb


# format sequence logo options
def formatOptions(logo, inOption, cs):
    if     inOption[2] : logo.color_scheme   = cs
    if not inOption[5] : logo.show_xaxis     = False
    if not inOption[6] : logo.show_yaxis     = False
    if not inOption[11]: logo.show_errorbars = False
    logo.unit_name          = inOption[0]
    logo.heat_scheme        = inOption[1]
    logo.show_colorkey      = inOption[3]
    logo.logo_title         = inOption[4]
    logo.yaxis_scale        = inOption[7]
    logo.xaxis_label        = inOption[8]
    logo.yaxis_label        = inOption[9]
    logo.yaxis_tic_interval = inOption[10]
    return logo
