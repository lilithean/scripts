#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
#
#              TITLE:           Calculate Lambda from Gamma and Omega
#             AUTHOR:           Xiaoyu Wang xwang224@buffalo.edu
#              USAGE:           clgo.py
#
#     REQUIRED FILES:
#       [fil_freq_gp]           produced by matdyn.x; default is freq.gp
#   [file_elph_gamma]           produced by matdyn.x; default is e.g. elph.gamma.6
#                               don't put X here as I'll add it later on
#    [fil_lambda_dat]           the file to read degauss and fermi level dos Nef
#
################################################################################
# Here you can set some parameters

# radius of lambda will times this number
magnitude = 1.
freq_tol =  15.

# input files
fil_freq_gp = 'freq.gp'
fil_elph_gamma = 'elph.gamma'
fil_lambda_dat = 'lambda.dat'

# constants
# 1 eV = 8065.5401069 cm-1
ev_wn = 8065.540107
wn_ev = 1./(ev_wn)
# 1 Ry = 13.605662 eV
ry_ev = 13.605662
# 1 eV = 241799.0504 GHz
ev_ghz = 241799.0504
ghz_ev = 1./ev_ghz
Pi = 3.1415927
e = 2.71828
wn_K = 1.4387863
wn_THz = 0.029979
################################################################################
# Below are main body of the program
################################################################################

import sys
from math import log

if len(sys.argv) >= 1:
    skip = True

def d_o_oqv(o, oqv, sigma=5.):
    return 1/((Pi*2)**0.5*sigma)*e**(-(o-oqv)**2/sigma**2)

# read lambda_dat
with open(fil_lambda_dat, 'r') as f:
    param = [[float(y) for y in x.split()] for x in f.readlines()[1:]]

# read fil_freq_gp
with open(fil_freq_gp, 'r') as f:
    Omega = [[float(y) for y in x.split()] for x in f.readlines()]

for i in range(2,10):
    # read fil_elph_gamma_X
    with open(fil_elph_gamma+'.'+str(i+1), 'r') as f:
        gdata = f.readlines()
    nbnd = int(gdata[0].split(',')[0].split()[2])
    nks = int(gdata[0].split()[4])
    nlk = int(nbnd/6 + bool(nbnd%6))
    Gamma = [
        [
            float(x) for x in ''.join(
                gdata[2+(nlk+1)*y:1+(nlk+1)*(y+1)]
            ).split()
        ] for y in range(nks)
    ]
    nef = param[i][4]/ry_ev
    lbda = param[i][1]

    # calculate lambda_qv
    Lambda = [
        [
            (Gamma[x][y]*ghz_ev)/(Omega[x][y+1]*wn_ev)**2/nef/Pi
            if Omega[x][y+1] > freq_tol else 0. for x in range(nks)
        ] for y in range(nbnd)
    ]
    #print(sum([sum([Lambda[y][x] for y in range(nbnd)]) for x in range(nks)]))

    # calculate a2F_qv
    a2F = [
        [
            (Gamma[x][y]*ghz_ev)/(Omega[x][y+1]*wn_ev)/nef/Pi/2.
            if Omega[x][y+1] else 0. for x in range(nks)
        ] for y in range(nbnd)
    ]

    # calculate A_qv
    A = [
        [
            (Gamma[x][y]*ghz_ev) * log(Omega[x][y+1]*wn_K)
            / lbda / nef / Pi / (Omega[x][y+1]*wn_ev)**2
            if Omega[x][y+1] > freq_tol else 0. for x in range(nks)
        ] for y in range(nbnd)
    ]

    # output bubble size files
    with open('a2F_band.'+str(i+1), 'w') as f:
        f.write(
            '\n\n'.join(
                '\n'.join(
                    "   %12.6f   %12.6f   %12.6f   " % (
                        Omega[x][0], Omega[x][y+1], a2F[y][x]
                    ) for x in range(nks)
                ) for y in range(nbnd)
            )
        )

    with open('gamma_band.'+str(i+1), 'w') as f:
        f.write(
            '\n\n'.join(
                '\n'.join(
                    "   %12.6f   %12.6f   %12.6f   " % (
                        Omega[x][0], Omega[x][y+1], Gamma[x][y]
                    ) for x in range(nks)
                ) for y in range(nbnd)
            )
        )

    with open('elph.lambda.'+str(i+1), 'w') as f:
        f.write(
            '\n\n'.join(
                '\n'.join(
                    "   %12.6f   %12.6f   %12.6f   " % (
                        Omega[x][0], Omega[x][y+1], Lambda[y][x]
                    ) for x in range(nks)
                ) for y in range(nbnd)
            )
        )

    with open('elph.lo.'+str(i+1), 'w') as f:
        f.write(
            '\n\n'.join(
                '\n'.join(
                    "   %12.6f   %12.6f   %12.6f   " % (
                        Omega[x][0], Omega[x][y+1],
                        Lambda[y][x] * Omega[x][y+1]
                    ) for x in range(nks)
                ) for y in range(nbnd)
            )
        )

    with open('elph.A.'+str(i+1), 'w') as f:
        f.write(
            '\n\n'.join(
                '\n'.join(
                    "   %12.6f   %12.6f   %12.6f   " % (
                        Omega[x][0], Omega[x][y+1], A[y][x]
                    ) for x in range(nks)
                ) for y in range(nbnd)
            )
        )

