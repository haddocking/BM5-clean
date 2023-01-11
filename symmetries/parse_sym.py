#!/usr/bin/env python3

import json
import sys
import os

with open(sys.argv[1]) as json_file:
    data = json.load(json_file)
    sym = []
    try:
        for s in data:
            sym.append(s['group'])
    except TypeError:
        pass

    print(','.join(sym), end='')

