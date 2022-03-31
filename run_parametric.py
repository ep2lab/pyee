'''
Creation Date: Tuesday April 27th 2021
Author: Pedro Jimenez, pejimene@ing.uc3m.es
##################################################
File: conv_study.py
Project: pyee
##################################################
Last Modified: Tuesday, April 27th 2021, 11:44:31 am
Modified By: Pedro Jimenez
##################################################
Code Description:
##################################################
Copyright (c) 2021 Plasma & Space Propulsion Team (EP2). Universidad Carlos III de Madrid.
'''
import run_parametric
import sys
import tools.parametric as parametric

def sim(case):

    n = [4,5,6,7,8,9,10]

    # This is an example you can use other functions or writte then in tools.parametric
    for i in n:
        parametric.nodes(case, int(1.25*2**i), int(1.25*2**(i-1)))

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    sim(*sys.argv[1:])