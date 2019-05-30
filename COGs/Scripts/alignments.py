### IMPORTS ###

import os
import re
from string import digits
import time
import numpy as np
from collections import defaultdict

### CHOOSE SOURCE FOLDER - Unhash the older of interest ###

# Combined
source_in = 'Unaligned'
source_out = 'Aligned'

### FUNCTIONS ###

def main():

    for root, dirs, filenames in os.walk(source_in):
        for f in filenames:

            #Define path for input/output files
            input = os.path.join(source_in, f)
            output = os.path.join(source_out, f)

            #Send input file to MAFFT, recieve aligned file
            cmd = "mafft-linsi %s > %s" % (input, output)
            os.system(cmd)

### RUN ###

if __name__ == '__main__':
    main()
