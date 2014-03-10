#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013 Junko Tsuji

# Generate and load dump data of logos.

import numpy as na
from corebio.matrix import Motif
from weblogolib import LogoData

# dump data to file
def dumpData(data, file):
    L = data.length
    f = open(file, "w")
    f.write("dim: " + str(data.length))
    f.write("," + str(len(data.alphabet)) + "\n")
    if data.entropy is not None:
        entropyStr = ",".join(map(str, data.entropy))
        f.write("ent: " + entropyStr + "\n")
    if data.weight is not None:
        weightStr = ",".join(map(str, data.weight))
        f.write("wgt: " + weightStr + "\n")

    for i in range(L):
        countStr = ",".join(map(str, data.counts.array[i]))
        f.write("cnt: " + countStr + "\n")
        if data.pvalue is not None:
            valueStr = ",".join(map(str, data.pvalue[i]))
            f.write("val: " + valueStr + "\n")
        if data.entropy_interval is not None:
            intervalStr = ",".join(map(str, data.entropy_interval[i]))
            f.write("ein: " + intervalStr + "\n")
        if data.odds_ratio is not None:
            oddsStr = ",".join(map(str, data.odds_ratio[i]))
            f.write("odd: " + oddsStr + "\n")
    if data.composition is not None:
        composStr = ",".join(map(str, data.composition))
        f.write("cmp: " + composStr + "\n")
    if data.max_value is not None:
        maxvalStr = str(data.max_value)
        f.write("max: " + maxvalStr + "\n")
    f.close()


# load dump data
def packData(file, alphabet, beg, end):
    length, chrlen = 0, 0
    counts = []
    entropy = []
    entropy_interval = []
    weight = []
    pvalue = []
    composition = []
    odds_ratio = []
    # read data from file
    for line in open(file, "r"):
        line = line.rstrip("\n")
        data = line.split()[1].split(",")
        if line.startswith("dim:"):
            length, chrlen = map(int, data)
        if line.startswith("ent:"):
            entropy = map(float, data)
        if line.startswith("wgt:"):
            weight = map(float, data)
        if line.startswith("cnt:"):
            counts.append(map(float, data))
        if line.startswith("val:"):
            pvalue.append(map(float, data[:chrlen]))
        if line.startswith("ein:"):
            entropy_interval.append(map(float, data))
        if line.startswith("odd:"):
            odds_ratio.append(map(float, data))
        if line.startswith("cmp:"):
            composition = map(float, data)
        if line.startswith("max:"):
            max_value = float(data[0])
    length = end - beg + 1
    if entropy == []:
        entropy = None
    else:
        entropy = entropy[beg-1:end]
    if entropy_interval == []:
        entropy_interval = None
    else:
        entropy_interval = entropy_interval[beg-1:end]
    if weight == []:
        weight = None
    else:
        weight = weight[beg-1:end]
    if pvalue == []:
        pvalue = None
    else:
        pvalue = pvalue[beg-1:end]
    if odds_ratio == []:
        odds_ratio = None
    else:
        odds_ratio = odds_ratio[beg-1:end]
    if max_value == []: max_value = None
    if composition == []: composition = None
    return LogoData(length, alphabet, Motif(array=counts[beg-1:end], alphabet=alphabet),
               entropy, entropy_interval, weight, pvalue, composition, odds_ratio, max_value)

