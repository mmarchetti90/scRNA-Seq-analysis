#!/usr/bin/env python3

"""
This script reformats the datasets/GSE115469_Data.csv.gz table.
"""

import pandas as pd

path = 'GSE115469_Data.csv.gz'

data = pd.read_csv(path, sep=',', skip_blank_lines=True)

data.columns = ['GeneSymbol'] + data.columns.to_list()[1:]

data.to_csv('datasets/GSE115469_Data.tsv.gz', compression='gzip', index=False, sep='\t')
