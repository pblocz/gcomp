#!/usr/bin/env python2
# coding=utf-8

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: [modified BSD](http://opensource.org/licenses/BSD-3-Clause)
'''

import json

import sympy as sp
from scipy.special import binom
from sympy.parsing.sympy_parser import parse_expr as parse


def gen_up_to_degree_and_save(deg, filename):
    x = sp.symbols('x')
    S = sp.symbols("a0:%d" % (deg + 1))

    res = {}
    for i in range(1, deg + 1):
        print i
        p = 0
        for k in range(i + 1):
            p += sp.binomial(i, k) * x**k * (1 - x)**(i - k) * S[k]

        res[i] = dict((k[0], v,) for k, v in p.as_dict().iteritems())

    def subs(expr):
        for i in range(deg + 1):
            expr = expr.subs(S[i], x**i)
            return expr

    coefs_l = []
    for i in range(1, deg + 1):
        print "degree", i
        c = [[0 for j in range(i + 1)] for k in range(i + 1)]
        for k, v in res[str(i)].iteritems():
            for deg, coef in subs(v).as_poly(x).as_dict().iteritems():
                deg, = deg
                c[int(k)][deg] = int(coef)
        coefs_l.append(c)

    with open(filename, "w") as f: json.dump(coefs_l, f)
