#!/usr/bin/env python3

import pandas as pd
import sys

df = pd.read_csv(str(sys.argv[1]), header=None)

minimum = df[1].min()

df[2] = (minimum / df[1])

headerList = ['ID', 'spikein_reads_mapped', 'scaling_factor']

df.to_csv('scaling_factors.csv', index=False, header=headerList)