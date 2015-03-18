# -*- encoding:utf-8 -*- 
from setuptools import setup
from distutils.core import  Extension

e_modules = [
    Extension('p3cbezier', sources = ['p3bezier/cbezier.c'], extra_compile_args=['-O2']),
]

setup (name = 'p3bezier',
       version = '1.0',
       url = 'https://github.com/pblocz/gcomp',
       author = "Pablo Cabeza & Diego González",
       author_email = "jcabeza@ucm.es | diegogon@ucm.es",
       description = 'This is a demo package',
       packages = ['p3bezier', 'p3bezier.binom'],
       
       data_files = [('p3bezier/binom/',['p3bezier/binom/binom.data'],)],
       ext_modules = e_modules,
       


       entry_points={
           'console_scripts': [
               'p3bezier = p3bezier.tester:main',
           ]
       },

)
