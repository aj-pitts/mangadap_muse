#!/usr/bin/env python3

#-----------------------------------------------------------------------
import sys
from os import remove
from mangadap.util.yanny import yanny
#-----------------------------------------------------------------------

def edit_dap_plan():

    if len(sys.argv) != 5:
        print('Usage: edit_dap_plan.py <plan file> <variable> <new value> <plan index>')
        raise Exception('Incorrect number of arguments!')

    print(sys.argv)

    par = yanny(filename=sys.argv[1])

    nplan = len(par['DAPPLAN'][sys.argv[2]])
    print('Number of plans read: {0}'.format(nplan))

    if sys.argv[4] == 'all':
        for i in range(0,len(par['DAPPLAN'][sys.argv[2]])):
            par['DAPPLAN'][sys.argv[2]][i] = sys.argv[3]
    else:
        if int(sys.argv[4]) < 0 or int(sys.argv[4]) >= nplan:
            raise Exception('{0} is not a valid index. Length is: {1}'.format(sys.argv[4], nplan))
        par['DAPPLAN'][sys.argv[2]][int(sys.argv[4])] = sys.argv[3]

    print(par['DAPPLAN'][sys.argv[2]])

    remove(sys.argv[1])
    par.write()

if __name__ == '__main__':
    edit_dap_plan()


