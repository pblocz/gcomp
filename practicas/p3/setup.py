from distutils.core import setup, Extension

e_modules = [
    Extension('cbezier', sources = ['bezier/cbezier.c']),
]

setup (name = 'bezier',
       version = '1.0',
       description = 'This is a demo package',
       packages = ['bezier'],
       ext_modules = e_modules,
)
