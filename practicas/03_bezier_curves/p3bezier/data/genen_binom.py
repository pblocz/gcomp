#!/usr/bin/env python2
# coding=utf-8

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: [modified BSD](http://opensource.org/licenses/BSD-3-Clause)
'''

from scipy.special import binom
import json
import os.path as path


def generate(N):
    l = []
    for n in range(N):
        for k in range(n + 1):
            l.append(binom(n, k))
    return l

with open(path.join(path.dirname(path.realpath(__file__)), "binom.data"), "r") as inputfile:
    table = json.load(inputfile)

if __name__ == "__main__":
    with open("binom.data", 'w') as outfile:
        json.dump(generate(500), outfile)
