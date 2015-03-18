#! /usr/bin/python -O
from scipy.special import binom
import json, os.path as path

def generate(N):
    l = []
    for n in range(N):
        for k in range(n+1):
            l.append(binom(n,k))
    return l

with open(path.join(path.dirname(path.realpath(__file__)),"binom.data"),"r") as inputfile:
    table = json.load(inputfile)

if __name__ == "__main__":
    with open("binom.data", 'w') as outfile:
        json.dump(generate(500), outfile)
