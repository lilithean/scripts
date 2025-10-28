#!/usr/bin/env python
#####################################################
#
# int_a2F.py [file:a2f.dos] <mu_star>
#
# 1. this script take QE output a2f.dos file which looks like:
# # omega total_a2f mode1_a2f mode2_a2f ...
# note that omega values given by QE are in Rydberg
#
# 2. this script output two files:
# - int_lambda.dat, which contains lambda(omega) values
# - int_omega.dat, which contains omega_ln(omega) values
#
#    Copyright(c) Xiaoyu Wang (xwang224@buffalo.edu)
#      Department of Wizardary and Alchemical Engineering
#      State University of New York at Buffalo, U.S.A.
#    Distributed under the terms of the MIT License.
#
#    https://ezwiki.wordpress.com/2021/12/19/electron-phonon-coupling-plot
#
#####################################################



import sys
#import argparse
from math import exp, log

RY2THZ = 3289.828                # Rydberg Constant to THz
THZ2K = 47.9924                 # THz to Kelvin
THZ2W = 33.35641                # THz to wavenumber
RY2W = RY2THZ * THZ2W           # Rydberg Constant to wavenumber
mustar = 0.1

def main():

    # read the data from a2F.dos.X file
    with open(sys.argv[1], 'r') as f:
        ls = f.readlines()
        dfreq = float(ls[-1].split()[5])
        freq = []
        a2f_tot = []
        a2f_mode = []
        for l in ls[5:-1]:
            freq.append(float(l.split()[0]))
            a2f_tot.append(float(l.split()[1]))
            a2f_mode.append([float(x) for x in l.split()[2:]])
        #for l in f[5:-1]:
        #    if len(l.split()[0]) == 8:
        #        freq.append(float(l.split()[0]))
        #        a2f_tot.append(float(l.split()[1]))
        #        a2f_mode.append([])
        #    else:
        #        a2f_mode[-1] += [float(x) for x in l.split()]


    if not dfreq:
        dfreq = freq[1] - freq[0]


    nfreq = len(freq) # number of frequency datapoints
    nmode = len(a2f_mode[0]) # number of modes

    gamma_tot = [-dfreq * (a2f_tot[0]*freq[0] + a2f_tot[-1]*freq[-1]) *RY2THZ**2.]
    for i in range(nfreq):
        gamma_tot.append(gamma_tot[-1]
                         + 2.0 * dfreq * a2f_tot[i]*freq[i] *RY2THZ**2.)

    ia2f_tot = [-dfreq * (a2f_tot[0] + a2f_tot[-1]) *RY2THZ]
    for i in range(nfreq):
        ia2f_tot.append(ia2f_tot[-1] + 2.0 * dfreq * a2f_tot[i] *RY2THZ)


    # lambda(omega) = 2 integral[0:omega]{domega a2F(omega)/omega}
    # for numerical integration, integral =
    #     sum[0:n-1]{domega/2 (a2f[i]/freq[i] + a2f[i+1]/freq[i+1])}
    #   = sum[0:n]{domega/2 a2f[i]/freq[i]} - domega/2 a2f[n]/freq[n]
    #   + sum[0:n]{domega/2 a2f[i]/freq[i]} - domega/2 a2f[0]/freq]0]
    #   = sum[0:n]{domega a2f[i]/freq[i]}        <----- this is the loop
    #   - domega/2 (a2f[n]/freq[n] + a2f[0]/freq[0])
    #       ^ this is the initial value of lambda_tot and lambda_mode

    lambda_tot = [-dfreq * (a2f_tot[0]/freq[0] + a2f_tot[-1]/freq[-1])]
    for i in range(nfreq):
        lambda_tot.append(lambda_tot[-1] + 2.0 * dfreq * a2f_tot[i]/freq[i])


    lambda_mode = [[-dfreq * (a2f_mode[0][x]/freq[0]
        + a2f_mode[-1][x]/freq[-1])
        for x in range(nmode)]]
    for i in range(nfreq):
        lambda_mode.append([lambda_mode[-1][x]
            + 2.0 * dfreq * a2f_mode[i][x]/freq[i]
            for x in range(nmode)])

    # omega_ln(omega) = exp(exponent(omega))
    # exponent(omega) =  2/lambda integral[0:omega]
    #                 {domega / omega a2f(omega) log(omega)}
    exp_tot = [-dfreq / lambda_tot[-1]
        * (a2f_tot[0] * log(abs(freq[0])*RY2THZ)/freq[0]
            + a2f_tot[-1] * log(freq[-1]*RY2THZ)/freq[-1])]
    for i in range(nfreq):
        exp_tot.append(exp_tot[-1] + 2./lambda_tot[-1]
            * a2f_tot[i] * dfreq * log(abs(freq[i])*RY2THZ)/freq[i])
    omega_ln_tot = [exp(x)*THZ2K for x in exp_tot]

    exp_mode = [[-dfreq / lambda_mode[-1][x]
        * (a2f_mode[0][x] * log(abs(freq[0])*RY2THZ)/freq[0]
            + a2f_mode[-1][x] * log(freq[-1]*RY2THZ)/freq[-1])
        for x in range(nmode)]]
    for i in range(nfreq):
        exp_mode.append([exp_mode[-1][x] + 2./lambda_mode[-1][x]
            * a2f_mode[i][x] * dfreq * log(abs(freq[i])*RY2THZ)/freq[i]
            for x in range(nmode)])
    omega_ln_mode = [[exp(y)*THZ2K for y in x] for x in exp_mode]


    # <omega_2> = 2/lambda_tot * int[0:omega]{domega a2f(omega) omega^(n-1)}
    omega_2 = [-dfreq / lambda_tot[-1] * RY2THZ**2.
        * (a2f_tot[0]*freq[0] + a2f_tot[-1]*freq[-1])]
    for i in range(nfreq):
        omega_2.append(omega_2[-1] + 2/lambda_tot[-1]
            * a2f_tot[i] * dfreq * freq[i]*RY2THZ**2.)

    omega_2 = omega_2[-1]**0.5*THZ2K

    L1 = 2.46 * (1. + 3.8*mustar)
    L2 = 1.82 * (1. + 6.3*mustar) * (omega_2/omega_ln_tot[-1])

    f1 = (1. + (lambda_tot[-1]/L1)**(1.5))**(1./3.)
    f2 = (1.
          + (omega_2/omega_ln_tot[-1] - 1.)
          * lambda_tot[-1]**2. / (lambda_tot[-1]**2. + L2**2.))

    Tc1 = (f1 * f2 * omega_ln_tot[-1] / 1.2
           * exp(-1.04*(1+lambda_tot[-1])
           / (lambda_tot[-1] - mustar - 0.62 * lambda_tot[-1] * mustar)))

    Tc2 = (omega_ln_tot[-1] / 1.2
           * exp(-1.04*(1+lambda_tot[-1])
           / (lambda_tot[-1] - mustar - 0.62 * lambda_tot[-1] * mustar)))

    print("lambda          = %12.4f" % lambda_tot[-1])
    print("omega_ln        = %12.4f" % omega_ln_tot[-1])
    print("omega_2         = %12.4f" % omega_2)
    print("f1              = %12.4f" % f1)
    print("f2              = %12.4f" % f2)
    print("Tc without f1f2 = %12.4f" % Tc2)
    print("Tc with f1f2    = %12.4f" % Tc1)
    print("lw              = %12.4f" % gamma_tot[-1])
    print("a2f             = %12.4f" % ia2f_tot[-1])

    with open("a2F.dat", "w") as f:
        for i in range(nfreq):
            f.write("%.6f       %.6f\n" % (freq[i]*RY2W, a2f_tot[i]))

    # output ALPHA2F.OUT file for eliash usage
    with open("ALPHA2F.OUT", 'w') as f:
        for i in range(nfreq):
            f.write("%.6f       %.6f\n" % (freq[i]/2, a2f_tot[i]))

    with open("int_lambda.dat", 'w') as f:
        for i in range(nfreq):
            f.write("  %10.6f   %10.6f   " % (freq[i]*RY2W, lambda_tot[i]))
            #f.write(" ".join(["%10.6f" % lambda_mode[i][x]
            #                    for x in range(nmode)]))
            f.write("\n")

    with open("int_a2f.dat", 'w') as f:
        for i in range(nfreq):
            f.write("  %.6f   %.6f  \n" % (freq[i]*RY2W, ia2f_tot[i]))

    with open("int_gamma.dat", 'w') as f:
        for i in range(nfreq):
            f.write("  %.6f   %.6f  \n" % (freq[i]*RY2W, gamma_tot[i]))

    with open("int_wln.dat", 'w') as f:
        for i in range(nfreq):
            f.write("  %.6f   %.6f  \n" % (freq[i]*RY2W, omega_ln_tot[i]))

if __name__ == "__main__":
    main()
