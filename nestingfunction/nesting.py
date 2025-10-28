#!/usr/bin/env python
import sys
import multiprocessing
from itertools import product
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator

class NestingCalculator(object):

    def __init__(self, fileigval, fildos):
        self.parallel = True
        self.nq1 = 20
        self.nq2 = 20
        self.nk1 = 160
        self.nk2 = 160
        self.etol = 0.01
        self.qtol = 0.2
        self.sigma_k = 0.04
        self.sigma_q = float(sys.argv[1])
        self.kmesh = [
            [  x,  y ]
            for x in np.linspace( .0, .5, self.nk1)
            for y in np.linspace( .0, .5, self.nk2)
        ]
        self.qmesh = [
            [  x,  y ]
            for x in np.linspace( .0, .5, self.nq1)
            for y in np.linspace( .0, .5, self.nq2)
        ]
        self.output_fil = 'nesting_func.dat'
        self.efbnds1 = []
        self.efbnds2 = []
        self.efkpts1 = []
        self.efkpts2 = []
        self.nestfunc = []
        self.peri = True
        #self.proc = 16

        with open(fileigval, 'r') as f:
            eigval = f.readlines()

        with open(fildos, 'r') as f:
            dos = f.readlines()

            self.efermi = float(dos[5].split()[3])
            self.nkpts = int(eigval[5].split()[1])
            self.nbnds = int(eigval[5].split()[2])

        self.kpts = [
            [float(x) for x in y.split()[:2]]
             for y in eigval[7::self.nbnds+2]
        ]
        self.bnds = []
        for i in range(self.nbnds):
            self.bnds.append(
                [float(
                    eigval[8+i+x*(self.nbnds+2)].split()[1]
                 )
                 - self.efermi
                for x in range(self.nkpts)]
            )

    def filt_bands(self):
        for i in range(self.nbnds):
            for j in range(len(self.kpts)):
                if abs(self.bnds[i][j]) < self.etol:
                    self.efbnds1.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts1.append(self.kpts[j])

                    self.efbnds1.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts1.append([
                        self.kpts[j][0], -self.kpts[j][1]
                    ])

                    self.efbnds1.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts1.append([
                        -self.kpts[j][0], self.kpts[j][1]
                    ])

                    self.efbnds1.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts1.append([
                        -self.kpts[j][0], -self.kpts[j][1]
                    ])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append(self.kpts[j])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append([
                         self.kpts[j][0],   -self.kpts[j][1]
                    ])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append([
                        -self.kpts[j][0],    self.kpts[j][1]
                    ])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append([
                        -self.kpts[j][0],   -self.kpts[j][1]
                    ])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append([
                        -self.kpts[j][0],  1-self.kpts[j][1]
                    ])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append([
                         self.kpts[j][0],  1-self.kpts[j][1]
                    ])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append([
                        1-self.kpts[j][0], 1-self.kpts[j][1]
                    ])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append([
                        1-self.kpts[j][0],   self.kpts[j][1]
                    ])

                    self.efbnds2.append(
                        Gauss(self.bnds[i][j], self.sigma_k)
                    )
                    self.efkpts2.append([
                        1-self.kpts[j][0],  -self.kpts[j][1]
                    ])



    def nest(self, q):
        s = 0.
        for i in range(len(self.efkpts1)):
            for j in range(len(self.efkpts2)):
                d_ij = np.array([
                    self.efkpts2[j][0] - self.efkpts1[i][0],
                    self.efkpts2[j][1] - self.efkpts1[i][1]
                ])
                d_q = np.linalg.norm(d_ij - q)
                if -self.qtol < d_q < self.qtol:
                    s += self.efbnds1[i] * self.efbnds2[j] \
                         * Gauss(d_q, self.sigma_q) \
                         / self.efdos1 / self.efdos2 \
                         / self.nk1**2. / self.nk2**2.
        return [s, q]

    def output_fs(self):
        with open("fromfs.dat", "w") as f:
            for i in self.efkpts1:
                f.write("%12.6f  %12.6f\n" % (i[0], i[1]))
        with open("tofs.dat", "w") as f:
            for i in self.efkpts2:
                f.write("%12.6f  %12.6f\n" % (i[0], i[1]))

    #def nest(self, i, j):
    #    q = [j[0]-i[0], j[1]-i[1]]
    #    return Gauss2d(self.qmesh, q, self.sigma_q)

    def calc_nest(self):
        if self.parallel:
            with multiprocessing.Pool() as pool:
                self.nestfunc = pool.map(self.nest, self.qmesh)

        else:
            self.nestfunc = [self.nest(x)
                for x in self.qmesh]

    def integral(self):
        s = 0.
        for i in self.nestfunc:
            s += i[0] * 0.5/self.nq1 * 0.5/self.nq2
        return s

    def output(self):
        with open(self.output_fil, 'w') as f:
            for i in self.nestfunc:
                #
                f.write("%12.6f %12.6f %18.9f\n" % (
                    i[1][0], i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    i[1][0], 1-i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    1-i[1][0], i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    1-i[1][0], 1-i[1][1], i[0]
                ))

                #
                f.write("%12.6f %12.6f %18.9f\n" % (
                    -i[1][0], i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    -i[1][0], 1-i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    -1+i[1][0], i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    -1+i[1][0], 1-i[1][1], i[0]
                ))

                #
                f.write("%12.6f %12.6f %18.9f\n" % (
                    i[1][0], -i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    1-i[1][0], -i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    i[1][0], -1+i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    1-i[1][0], -1+i[1][1], i[0]
                ))

                #
                f.write("%12.6f %12.6f %18.9f\n" % (
                    -i[1][0], -i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    -1+i[1][0], -i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    -i[1][0], -1+i[1][1], i[0]
                ))
                f.write("%12.6f %12.6f %18.9f\n" % (
                    -1+i[1][0], -1+i[1][1], i[0]
                ))


    def interpolate(self):
        with open("bnd.dat", 'w') as f:
            for i in range(len(self.kpts)):
                f.write("%12.6f %12.6f %19.8f\n" % (
                     self.kpts[i][0],
                     self.kpts[i][1],
                     self.bnds[0][i]))
        intp_bnds = []
        for x in self.bnds:
            interp = CloughTocher2DInterpolator(
                self.kpts, x
            )
            intp_bnds.append(interp(self.kmesh))

        self.kpts = self.kmesh
        self.bnds = intp_bnds

        with open("interp_bnd.dat", 'w') as f:
            for i in range(len(self.kmesh)):
                f.write("%12.6f %12.6f %19.8f\n" % (
                     self.kmesh[i][0],
                     self.kmesh[i][1],
                     intp_bnds[0][i]))

    def dosef(self):
        self.efdos1 = sum(self.efbnds1) \
               * (0.5/self.nk1) * (0.5/self.nk2)
        self.efdos2 = sum(self.efbnds2) \
               * (0.5/self.nk1) * (0.5/self.nk2)

    def check_mesh(self):
        return all(x in self.qmesh for x in self.kmesh)

    def run(self):
        print('number of cores', multiprocessing.cpu_count())
        print('number of k-points', self.nkpts)
        print('number of bands', self.nbnds)
        print('fermi energy', self.efermi)
        print('energy tolerence within Ef', self.etol)
        print('interpolate eigenvalues...')
        self.interpolate()
        print('filter eigenvalues near Ef')
        self.filt_bands()
        print('number of points near Ef', len(self.efkpts1))
        self.output_fs()
        self.dosef()
        print('DOS at EF (from)', self.efdos1)
        print('DOS at EF (to)', self.efdos2)
        print('number of q-vecs to estimate', len(self.qmesh))
        print('calculate nesting function...')
        self.calc_nest()
        print('nesting calculation done')
        print('integral of nesting func is',
               self.integral())
        print('gamma point nesting func', self.nestfunc[1][0])

        self.output()

def Gauss(x, sigma, x0=0.):
    return np.exp(
        -((x-x0)**2 / (2.*sigma**2))
    ) / sigma / np.sqrt(2. * np.pi)

def Gauss2d(xy, xy0, sigma=0.01):
    return np.array([
        [  (np.exp(-((x[0]-xy0[0])**2 / (2.*sigma**2)))
          + np.exp(-((x[1]-xy0[1])**2 / (2.*sigma**2))))
         / sigma**2. / 2. / np.pi
         for x in xy]
    ])

def main():
    Nest = NestingCalculator('EIGENVAL', 'DOSCAR')
    Nest.run()

if __name__ == '__main__':
    main()
