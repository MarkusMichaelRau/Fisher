import numpy as np

import sys

if (__name__ == '__main__'):

    args = sys.argv[1:]

    if len(args) != 1:
        print 'Invalid number of command line arguments in pad_dndz.py'
        sys.exit(0)

    data = np.loadtxt(args[0])
    data[0, 1:] = 0.0
    data[-1, 1:] = 0.0
    np.savetxt(X=data, fname=args[0])
