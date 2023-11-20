# Maxent.py
# @author Yanfei Tang
# modified by Thomas A. Maier

import sys
import os
import numpy as np
import time
import _pickle as pickle
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.optimize import minimize
from numpy.linalg import inv

def gaussian(w, mu = 0.0, sig = 4.0):
    """
    A gaussian distribution.
    """
    return 1.0/sig/np.sqrt(2* np.pi) * np.exp(-np.power(w - mu, 2.)/ 2. / np.power(sig, 2.))

def straightline(w):
    """
    A straight line.
    """
    return np.ones(len(w))/(w[-1] - w[0])

class Maxent(object):
    """
    Maximum entropy method of analytical continuation of Imaginary frequency Green's function data.

    Attributes:
    ----------
    kernel:
        fermionic_frequency_withRealPart for fermionic frequency kernel if real part of data != 0
        fermionic_frequency_noRealPart   for fermionic frequency kernel if real part of data  = 0
        fermionic_time                   for fermionic imag. time kernel
        bosonic_frequency                for bosonic frequency kernel
        bosonic_time                     for bosonic imag. time kernel
    wn: ndarray
        container for eith Matsubara frequency, e.g. fermionic case wn = (2n+1)/pi/beta,
        or time
    G: ndarray, dtype = 'complex'
        imaginary frequency Green's function.
    K: ndarray, dtype = 'complex'
        Kernel.
    w: ndarray
        real frequency.
    aveG: ndarray, dtype = 'complex'
        average of imaginary frequency Green's function.
    stdG: ndarray, dtype = 'complex'
        standard deviation of imaginary frequency Green's function.
    numMfre: int
        number of Matsubara frequency or time steps, set from the input data
    numbins: int
        number of bins, e.g. 20 sets of QMC data.
    numRfre: int
        number of real frequency.
    specF: ndarray
        spectral function.
    defaultM: ndarray
        default model.
    wmin: float
        minimum real frequency.
    wmax: float
        maximum real frequency.
    alpha: float
        tuning parameter to adjust the trade-off between the fit of the data and the nearness of the default model.
    dw: float
        the discrete space of w.
    tol: float
        tolerence of the minimization.
    std: tuple (bool, float)
        whether considering the standard deviation into the chi sqaure computation.
    prob: float
        the probability for this alpha.
    alphamin and alphamax and dalpha: float
        the period of alpha, and the interval of alpha
    alphas: ndarray
        the alpha that are examed.
    allSpecFs: ndarray
        all the spectral functions for different alpha
    allProbs: ndarray
        probabilities for each different alpha
    aveSpecFs: ndarry
        the average of the spectral functions
    mu: float
        parameter in Levenberg-Marquart algorithms.
    maxIteration: int
        Maximum iteration in Lenvenberg-Marquart algorithm.

    """


    # def __init__(self, filename = "./Gr_7.txt", beta=10, kernel = "fermionic_frequency_withRealPart", numRfre = 501, wmin = -15, wmax = 15, defaultModel = "gaussian", dCenter = 0, gaussianSigma=4, tol = 1e-7, std = (True, 1.0), alphamin = -2, alphamax = 1, numAlpha = 100, draw = True):
    def __init__(self, filename, beta, kernel, numRfre, wmin, wmax, defaultModel, dCenter, gaussianSigma, tol, std, alphamin, alphamax, numAlpha, draw):
    
        """
        Contructor
        """
        print("wmin %s" %wmin)
        print("beta %s" %beta)
        self.start_time = time.time()
        self.numRfre = numRfre
        self.wmin = wmin
        self.wmax = wmax
        self.kernel = kernel
        self.beta=beta
        self.tol = tol
        self.std = std
        self.alphamin = alphamin
        self.alphamax = alphamax
        self.numAlpha = numAlpha
        self.draw = draw
        self.mu = 0.1
        self.maxIteration = 1e4
        self.allSpecFs, self.allProbs = [], []

        self.readfile(filename)
        # self.jackknife()

        self.calMeanAndStd()
        self.createKernel()
        self.rotation()
        self.compUXiV()
        if defaultModel == 'gaussian':
            self.specF = gaussian(self.w,mu=dCenter,sig=gaussianSigma)/self.normalize(gaussian(self.w,mu=dCenter,sig=gaussianSigma))
            self.defaultM = self.specF
        elif defaultModel == 'doubleGaussian':
            m1 = gaussian(self.w,mu=-dCenter,sig=gaussianSigma)/self.normalize(gaussian(self.w,mu=-dCenter-8,sig=gaussianSigma))
            m2 = gaussian(self.w,mu=dCenter,sig=gaussianSigma)/self.normalize(gaussian(self.w,mu=dCenter+8,sig=gaussianSigma))
            self.specF = 0.5*(m1 + m2)
            self.defaultM = self.specF
        elif defaultModel == 'straightline':
            self.specF, self.defaultM = straightline(self.w), straightline(self.w)
        elif defaultModel == 'input':
            # data = np.loadtxt(sys.argv[2])
            data = np.loadtxt("defaultModel.txt")
            if data.shape[0] != self.numRfre:
                sys.exit("Number of frequencies in default model does not match numRfre; exiting!")
            self.specF, self.defaultM = data[:, 1], data[:, 1]
        else:
            print("Usage: Maxent(filename, numRfre, wmin, wmax, defaultModel, tol, std, alphamin, alphamax, numAlpha)")
            print("defaultModel must 'g' or 's'!")
            sys.exit(0)
        print("Use --- %s seconds --- before optimization..." %(time.time() - self.start_time))

    def getAllSpecFs(self):
        """
        Compute all the spectral functions by looping all the alphas.
        """
        self.alphas = np.logspace(self.alphamin, self.alphamax, self.numAlpha)
        # Use uniformed space integral for \int P(alpha)dalpha;
        self.dalpha = self.alphas[1:] - self.alphas[:-1]
        self.dalpha /= 2.0
        self.dalpha = np.insert(self.dalpha, 0, 0.0) + np.append(self.dalpha, 0.0)
        for self.alpha in self.alphas:
            self.getSpecF()
            self.calProb()
            self.allSpecFs.append(self.specF)
            self.allProbs.append(self.prob)
            print("Finish alpha = %s.\n" %(self.alpha))

        self.allProbs = np.array(self.allProbs)
        self.allSpecFs = np.array(self.allSpecFs)
        self.allProbs = self.allProbs/np.sum(self.allProbs*self.dalpha)
        self.aveSpecFs = np.dot(self.allSpecFs.transpose() * self.dalpha, self.allProbs)
        print("Optimization ends. Use --- %s seconds ---" %(time.time() - self.start_time))

    def saveObj(self,dirname):
        """
        save the object as .p file
        And plot the figure.
        rFre:
            real frequency
        tSpec:
            true spectral function

        """
        
        # if sys.argv[1][-1] == "G":
        #     dirname = sys.argv[1][:-1]
        # else:
        #     dirname = sys.argv[1]
        #---------save object pickle-----------------
        #ofile = open(dirname + "object.p", 'w')
        #pickle.dump(self, ofile, -1)
        #ofile.close()
        #---------save dict pickle-------------------
        dict_data = {"wn": self.wn,
                     "w": self.w,
                     "aveG": self.aveG,
                     "stdG": self.stdG,
                     "numMfre": self.numMfre,
                     "numbins": self.numbins,
                     "numRfre": self.numRfre,
                     "defaultM": self.defaultM,
                     "alphas": self.alphas,
                     "allSpecFs": self.allSpecFs,
                     "allProbs": self.allProbs,
                     "aveSpecFs": self.aveSpecFs,
        }

        # if os.path.isfile("./" + dirname + "A"):
        #     rFre, tSpecF = [], []
        #     with open(dirname+ "A", "r") as ifile:
        #         for line in ifile:
        #             term = line.split()
        #             rFre.append(float(term[0]))
        #             tSpecF.append(float(term[1]))
        #     rFre, tSpecF = np.array(rFre), np.array(tSpecF)
        #     dict_data["rFre"] = rFre
        #     dict_data["tSpecF"] = tSpecF
        #     plt.plot(rFre, tSpecF, "r--", alpha = 0.5, label = "TestFunc")
        #     ifile.close()

        # ofile = open(dirname + ".p", 'wb')
        # pickle.dump(dict_data, ofile, -1)
        # ofile.close()

        result = open(dirname + "maxent.dat", 'w')
        for i in range(len(self.w)):
            result.write(str(self.w[i]) + '\t' + str(self.aveSpecFs[i]) + "\n")
        result.close()

        # result = open(dirname + "Palpha.dat", 'w')
        # for i in range(len(self.alphas)):
        #     result.write(str(self.alphas[i]) + '\t' + str(self.allProbs[i]) + "\n")
        # result.close()

        # if self.draw:
        #     plt.plot(self.w, self.aveSpecFs,  alpha = 0.8, label = "Maxent")
        #     plt.xlabel(r"$\omega$")
        #     plt.ylabel(r"$A(\omega)$")
        #     plt.legend()
        #     plt.savefig("./" + dirname + "Comparison.pdf")
        #     plt.show()


    def readfile(self, filename):
        print("Reading file from ", filename)
        """
        read data from file:
        file format:
        #   1   2   3   4   5  ... column
        1  wn ReG ImG ReG ImG ...  ImG
        2  wn ... ... ... ... ...  ImG
        ...
        row

        Note:   Reads only positive Mats. frequencies and either both real and imaginary parts of data when
                self.kernel = 'fermionic_frequency_withRealPart' or else only single column
        """

        data = np.loadtxt(filename)
        self.wn = data[:,0]
        self.numMfre = self.wn.shape[0]

        step = 1
        if self.kernel in ('fermionic_frequency_withRealPart','fermionic_frequency_noRealPart','bosonic_frequency_noRealPart','bosonic_frequency_withRealPart') : 
            step = 2
            
        if step==1:
            self.G = data[:,1:]
        else: # Data contains real and imaginary part
            ReG = data[:,1::2]
            ImG = data[:,2::2]
            if self.kernel in ('fermionic_frequency_noRealPart','bosonic_frequency_noRealPart'):
                self.G = ImG
            else:
                numSamples = ReG.shape[1]
                self.G = np.zeros((2*self.numMfre,numSamples))
                self.G[0::2,:] = ReG[:,:]
                self.G[1::2,:] = ImG[:,:]
        self.numGfre = self.G.shape[0]
        print("file ",filename," contains ",self.numMfre, " positive frequencies/time slices and ", self.G.shape[1] ," samples.")
        if self.G.shape[1] < 2*self.numMfre:
            print("Shut Down immediately! number of samples must be larger than or equal to 2 times of number of Matsubara frequencies.")
            sys.exit(0)


    def jackknife(self):
        """Jackknife procedure"""
        self.numMfre, self.numbins = self.G.shape
        idx = np.arange(self.numbins)
        est = np.zeros((self.numMfre,self.numbins),dtype=type(self.G[0][0]))
        for i in range(self.numbins):
            est[:,i] = np.sum(self.G[:,idx!=i],axis=1)/float(self.numbins-1)
        self.G = est

    def calMeanAndStd(self):
        """
        From data (G) to calculate the average, standard deviation and covariance matrix.
        """
        self.numMfre, self.numGfre,self.numbins = self.wn.shape[0], self.G.shape[0],self.G.shape[1]
        self.aveG = self.G.mean(1)
        A = np.array([[self.aveG[i]]*self.numbins for i in range(self.numGfre)]).reshape(self.numGfre, self.numbins) - self.G
        self.stdG = np.sqrt(np.sum(A * A, axis = 1)/(self.numbins - 1))
        self.cov = np.zeros([self.numGfre, self.numGfre], dtype = 'float64')
        for l in range(self.numGfre):
            for k in range(self.numGfre):
                a = 0.0
                for j in range(self.numbins):
                    a += (self.aveG[l] - self.G[l][j] ) * (self.aveG[k] - self.G[k][j])
                self.cov[l][k] = a/(self.numbins-1)


    def createKernel(self):
        self.numMfre = len(self.wn)
        self.numGfre = self.G.shape[0]
        self.w = np.linspace(self.wmin, self.wmax, self.numRfre)
        self.dw = self.w[1:] - self.w[:-1]
        self.dw = np.append(self.dw, self.dw[0])
        self.dw[0] /= 2.0
        self.dw[-1] /= 2.0
        self.K = np.zeros([self.numGfre, self.numRfre], dtype = 'float64')
        if self.kernel == "fermionic_frequency_withRealPart":
            """
            The fermionic kernel is:
            #         1
            # K = ---------.
            #     i*wn - w
            """
            for n in range(self.numMfre):
                for m in range(self.numRfre):
                    c1 = 1.0/(-self.w[m] + self.wn[n] * 1j)
                    self.K[2*n][m]   = np.real(c1)
                    self.K[2*n+1][m] = np.imag(c1)

        if self.kernel == "fermionic_frequency_noRealPart":
            """
            The fermionic kernel is:
            #         1
            # K = ---------.
            #     i*wn - w
            """
            for n in range(self.numMfre):
                for m in range(self.numRfre):
                    c1 = 1.0/(-self.w[m] + self.wn[n] * 1j)
                    self.K[n][m] = np.imag(c1) # Data and kernel contain only imaginary parts

        if self.kernel == "fermionic_time":
            """
            The fermionic kernel is:
            #         exp(-tau*w)
            # K = - ------------------
            #        1 + exp(-w*beta)
            """
            beta = self.beta
            for n in range(self.numMfre):
                for m in range(self.numRfre):
                    w = self.w[m]
                    tau = self.wn[n]
                    self.K[n][m] = np.exp(-w*tau)/(1.0+np.exp(-w*beta))

        if self.kernel == "bosonic_frequency_noImagPart":  # When the imag part = 0, this means the spectrum is ph-symmetric.
                """
                #         w^2
                # K = ------------
                #     wn^2 + w^2
                """
                for n in range(self.numMfre):
                    for m in range(self.numRfre):
                        if (self.wn[n]==0) & (self.w[m]==0):
                            self.K[n][m] = 1.0
                        else:
                            self.K[n][m] = self.w[m]**2/(self.w[m]**2 + self.wn[n]**2)

        if self.kernel == "bosonic_frequency_withRealPart":
                """
                #         w
                # K = ------------
                #     iwm + w
                """
                for n in range(self.numMfre):
                    for m in range(self.numRfre):
                        if (self.wn[n]==0) & (self.w[m]==0):
                            self.K[2*n][m]   = 1.0
                            self.K[2*n+1][m] = 0.0
                        else:
                            c1 = self.w[m]/(1j*self.wn[n] + self.w[m])
                            self.K[2*n][m]   = np.real(c1)
                            self.K[2*n+1][m] = np.imag(c1)

        if self.kernel == "bosonic_time":
            """
            # If data is not symmetric, i.e. G(tau) != G(beta-tau), i.e. G(iwm) != -G(-iwm)
            # 
            #         w * exp(-tau*w)
            # K =   ------------------
            #        1 - exp(-w*beta)
            """
            beta = self.beta
            for n in range(self.numMfre):
                for m in range(self.numRfre):
                    w = self.w[m]
                    if w==0.0:
                        self.K[n][m] = 1./beta
                    else:
                        tau = self.wn[n]
                        self.K[n][m] = w*np.exp(-w*tau)/(1.0-np.exp(-w*beta))

        if self.kernel == "bosonic_time_symmtric":
            """
            # If data is symmetric, i.e. G(tau) = G(beta-tau), i.e. G(iwm) = -G(-iwm)
            # 
            #             exp(-w*tau)+exp(-w(beta-tau))
            # K = 1/2*w  --------------------------------
            #                    1 - exp(-w*beta)
            """
            beta = self.beta
            for n in range(self.numMfre):
                for m in range(self.numRfre):
                    w = self.w[m]
                    if w==0.0:
                        self.K[n][m] = 1./beta
                    else:
                        tau = self.wn[n]
                        self.K[n][m] = 0.5*w*(np.exp(-w*tau)+np.exp(-w*(beta-tau)))/(1.0-np.exp(-w*beta))


    def rotation(self):
        """
        Rotate the kernel and data into diagonal representation.
        """
        w, v = np.linalg.eigh(self.cov)
        vt = v.transpose().conjugate()
        self.stdG = np.sqrt(w)
        self.K = np.dot(vt, self.K)
        self.aveG = np.dot(vt, self.aveG)


    def compUXiV(self):
        """
        Do the sigular value decomposition to the kernel matrix K.
        """

        self.U, self.Xi, self.Vt = np.linalg.svd(self.K.transpose(), full_matrices = 0)
        self.rank = np.linalg.matrix_rank(self.K.transpose())
        self.Xi = np.diag(self.Xi[:self.rank])
        self.Vt = np.real(self.Vt[:self.rank, :])
        self.V = self.Vt.transpose()
        self.U = np.real(self.U[:, :self.rank])
        self.Ut = self.U.transpose()


        self.M = np.dot(np.dot(self.Xi, self.Vt), np.diag(1.0/self.stdG/self.stdG))
        self.M = np.dot(self.M, self.V)
        self.M = np.dot(self.M, self.Xi)

    def chiSquare(self, specF):
        """
                   1
        \chi^2 = ----- |(aveG - K * A)/\sigma|^2
                   N
        """
        delta = self.aveG - np.dot(self.K, specF * self.dw)
        return np.real( np.sum( delta * delta/self.stdG/self.stdG )/self.numGfre )

    def restoreG(self, specF):
        """
        From the spectral function to find out the tilted G: G = K * A
        """
        return np.dot(self.K, specF * self.dw)

    def normalize(self, specF):
        """
        Find out the normalization facotr \int A(w) * dw
        """
        return np.sum( specF * self.dw )

    def objective(self, specF):
        """
        Q = 1/2 \chi^2 - \alpha * S
        considering the standard deviation or not.
        """
        delta = self.aveG - np.dot(self.K, specF * self.dw)
        if self.std[0]:
            return np.real(np.sum( delta*delta/self.stdG/self.stdG ))/2.0 + \
                self.alpha * np.sum((specF * np.log(np.abs((specF)/self.defaultM)) - specF +  self.defaultM) * self.dw)
        else:
            return np.real(np.sum( delta*delta/self.std[1]/self.std[1] ))/2.0 + \
                self.alpha * np.sum((specF * np.log(np.abs((specF)/(self.defaultM))) - specF +  self.defaultM) * self.dw)

    def getSpecF(self):
        """
        using Bryan's method (R. K. Bryan, Eur. Biophys. J., 18 (1990) 165) to minimize the objective function to get the spectral function depending on alpha.
        """
        iteration = 0
        n = min(self.numMfre, self.numRfre)
        btemp = np.zeros(self.rank)
        self.specF = self.defaultM
        Qold = self.objective(self.defaultM)

        while True:
            iteration += 1
            T = np.dot(np.dot(self.Ut, np.diag(self.specF)), self.U)
            deri = -1.0/self.stdG/self.stdG * (self.aveG - self.restoreG(self.specF))
            g = np.dot(np.dot(self.Xi, self.Vt), deri)
            LHS = (self.alpha + self.mu) * np.diag(np.ones(self.rank)) + np.dot(self.M, T)
            RHS = -self.alpha * btemp - g
            deltab = np.dot(inv(LHS), RHS)
            criteria = np.dot(deltab, np.dot(T, deltab))

            if criteria < 0.2*sum(self.defaultM):
                if iteration > self.maxIteration:
                    print("Exceeds maximum iteration in Levenberg-Marquart algorithms, exits. Increase tolerance.")
                    break

                btemp += deltab
                al = np.dot(self.U, btemp)
                self.specF = self.defaultM * np.exp(al)
                Qnew = self.objective(self.specF)
                if abs(Qnew - Qold)/Qold < self.tol:
                    print("{0} iterations in Levenberg-Marquart algorithms. Function evaluted: {1}, it exits properly.".format(iteration, Qnew))
                    break
                Qold = Qnew
                continue

            else:
                self.mu *= 1.3
                self.specF = self.defaultM
                Qold = self.objective(self.defaultM)
                btemp = np.zeros(self.rank)
                print("parameter \mu is too small in the Levenberg-Marquart algorithms.")
                print("\mu is now adjusted to %s" %(self.mu))



    def calProb(self):
        """
        Compute the probablity for this paticular alpha. This probablity is not normalized.
        """
        if self.std[0]:
            cov = np.diag(self.stdG * self.stdG)
        else:
            cov = np.diag(np.ones(self.numGfre) * np.ones(self.numGfre) * self.std[1] * self.std[1])
        mat_a = np.dot(self.K.transpose(), inv(cov))
        mat_a = np.dot(mat_a, self.K)
        vec_a = np.sqrt(np.abs(self.specF))
        imax, jmax = mat_a.shape
        mat_b = np.zeros((imax, jmax)) + 1j * 0
        for i in range(0, imax):
            for j in range(0, jmax):
                mat_b[i][j] = vec_a[i] * mat_a[i][j] * vec_a[j]
        S = np.linalg.eigvalsh(mat_b)
        expo = np.exp(-self.objective(self.specF))
        prod = np.prod(self.alpha/(self.alpha+S))


        self.prob = np.sqrt( prod ) * expo/self.alpha

        if np.isnan(self.prob):
            self.prob = 0.0




if __name__ == "__main__":
    """
    filename: the file that stores imaginary-frenquency Green's function data
    numMfre: number of Matsubara frequency used. (this number should be less than the rows in the file)
    numRfre: the number of grid for spectral function A(w).
    wmin: minimum real frequency
    wmax: maximun real frequency
    defaultModel: this parameter can be 'gaussian', 'straightline' or 'input'.
    tol: tolerance for minimization. 1e-12 for SLSQP; 1e-5 for Bryan
    std: whether or not use the standard deviation.
    alphamin: value of minimun alpha in log space.
    alphamax: value of maximum alpha in log space
    numAlpha: number of alphas.
    minimizer: "SLSQP" or "Bryan".
    draw: whether or not draw the Maxent result graph.
    """
    # Model = Maxent(filename = sys.argv[1], statistics='B', numMfre = 48, numRfre = 100, wmin = -20, wmax = 20, defaultModel = 'gaussian', tol = 1e-3, std = (True, 1.0), alphamin = -1, alphamax = 2, numAlpha = 100, minimizer = "Bryan", draw = True)
    Model = Maxent(filename = sys.argv[1], beta=10., kernel = "bosonic_time_symmetric", numRfre = 201, wmin = -10, wmax = 10, defaultModel = "gaussian", dCenter=0, gaussianSigma=4, tol = 1e-7, std = (True, 1.0), alphamin = -2, alphamax = 1, numAlpha = 100, draw = True)
    Model.getAllSpecFs()
    Model.saveObj(sys.argv[2])





