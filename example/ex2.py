#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Nov  3 16:29:34 2017

Zewei Song
BGI-Shenzhen
songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
from rrna_parser import data_io
import importlib

importlib.reload(data_io)
# from rrna_parser import data_io/This line is not necessary
a = 'ACACTG'
b = data_io.revcomp(a)
print(a,b)


from __future__ import print_function
from __future__ import division
import rrna_parser.data_io as data_io
import importlib

importlib.reload(data_io)
#import rrna_parser.data_io as data_io
a = 'ACACTG'
b = data_io.revcomp(a)
print(a,b)