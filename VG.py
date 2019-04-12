from JuZhong import JuZhongPrice
from math import log, sqrt, exp, pi, gamma
from numpy import sign
from scipy.integrate import quad
from scipy.stats import norm
from scipy.special import kv
import numpy as np
import matplotlib.pyplot as plt


def VGEuro(S, K, nu, theta, sigma, r, q, T, phi):
    '''
    phi = 1: call; phi = -1: Put
    '''

    c = T/nu
    omega = log(1 - nu * theta - (sigma**2) * nu / 2)/nu
    zeta = (log(S/K) + omega * T)/sigma
    v = 1 - nu * (theta + sigma**2 / 2)
    a1 = zeta * sqrt(v/nu)
    b1 = (theta + sigma**2)*sqrt(nu/v) / sigma
    a2 = zeta * sqrt(1/nu)
    b2 = theta * sqrt(nu)/sigma
    euroCall = S * exp(-q * T) * Phi(a1, b1, c) - K * exp(-r * T) * Phi(a2, b2, c)


    if phi == 1:
        return euroCall
    else:
        return euroCall - (S * exp(-q*T) - K * exp(-r*T))

def Phi(a, b, gammaPhi):
    u = b / sqrt(2+b**2)
    c = abs(a) * sqrt(2 + b**2)
    phi1 = hyperGM(gammaPhi, 1 - gammaPhi, 1 + gammaPhi, (1+u)/2, -sign(a)*c*(1+u))
    phi2 = hyperGM(1+gammaPhi, 1 - gammaPhi, 2+gammaPhi, (1+u)/2, -sign(a)*c*(1+u))
    param1 = (c**(gammaPhi+0.5) * exp(sign(a)*c) * (1+u)**gammaPhi) / (sqrt(2*pi)*gamma(gammaPhi)*gammaPhi)
    param2 = (c**(gammaPhi+0.5) * exp(sign(a)*c) * (1+u)**(1+gammaPhi)) / (sqrt(2*pi)*gamma(gammaPhi)*(gammaPhi+1))
    return param1 * phi1 * kv(gammaPhi+0.5, c) - sign(a) * param2 * kv(gammaPhi-0.5, c) * phi2 + sign(a)*param1*kv(gammaPhi-0.5, c)*phi1

def hyperGM(alpha, beta, gam, x, y):
    intSol = 0
    integrandL = 0
    for i in range(1, 100001):
        deltaU = 1 / 100000
        u = deltaU * i
        if i != 1:
            integrandL = integrandH
        if 1 - u != 0:
            integrandH = (u**(alpha-1)) * ((1-u)**(gam - alpha - 1)) * ((1 - u * x)**(-beta)) * exp(u*y)
        else:
            integrandH = (u ** (alpha - 1)) * ((1 + deltaU/2 - u) ** (gam - alpha - 1)) * ((1 - u * x) ** (-beta)) * exp(u * y)
        intSol += (integrandL + integrandH) * deltaU / 2
    return (intSol * gamma(gam)) / (gamma(alpha) * gamma(gam - alpha))

def VGPrice(S, K, nu, theta, sigma, r, q, T, phi):
    def expint(lowerBound):
        def integrand(y):
            return exp(-y) / y
        return quad(integrand, lowerBound, np.inf)[0]
    lambdaP = sqrt((theta**2)/(sigma**4) + 2/(sigma**2 * nu)) - theta/(sigma**2)
    lambdaN = lambdaP + 2 * theta/(sigma**2)

    epsilon = S/8000
    sigmaEpsilon = sqrt((1 - exp(-lambdaP*epsilon) * (epsilon*lambdaP + 1)) / (nu * (lambdaP**2)) + (1 - exp(-lambdaN*epsilon) * (epsilon*lambdaN + 1)) / (nu * (lambdaN**2)))
    omegaEpsilon = expint(epsilon * lambdaP) / nu - expint(epsilon * (lambdaP - 1)) / nu + expint(epsilon * lambdaN) / nu - expint(epsilon * (lambdaN + 1)) / nu
    return JuZhongPrice(S, K, r, q-omegaEpsilon, sigmaEpsilon, T, phi)[-1], VGEuro(S, K, nu, theta, sigma, r, q, T, phi)
    #return JuZhongPrice(S, K, r, q-omegaEpsilon, sigmaEpsilon, T, phi)[0]



def paramEva(S,nu, theta, sigma, r, q, T):
    def expint(lowerBound):
        def integrand(y):
            return exp(-y) / y
        return quad(integrand, lowerBound, np.inf)[0]
    lambdaP = sqrt((theta**2)/(sigma**4) + 2/(sigma**2 * nu)) - theta/(sigma**2)
    lambdaN = lambdaP + 2 * theta/(sigma**2)
    epsilonList = [i/100 for i in range(200)]
    sigmaList = []
    omegaList = []
    for epsilon in epsilonList:
        sigmaEpsilon = sqrt((1 - exp(-lambdaP*epsilon) * (epsilon*lambdaP + 1)) / (nu * (lambdaP**2)) + (1 - exp(-lambdaN*epsilon) * (epsilon*lambdaN + 1)) / (nu * (lambdaN**2)))
        omegaEpsilon = expint(epsilon * lambdaP) / nu - expint(epsilon * (lambdaP - 1)) / nu + expint(epsilon * lambdaN) / nu - expint(epsilon * (lambdaN + 1)) / nu
        sigmaList.append(sigmaEpsilon)
        omegaList.append(omegaEpsilon)
    #plt.plot(epsilonList, sigmaList, 'o')
    plt.plot(epsilonList, omegaList, 'or')
    plt.title('T = {}'.format(T))
    plt.xlabel('epsilon')
    #plt.ylabel('sigma')
    plt.ylabel('omega')
    plt.legend(loc='upper left')
    plt.show()

#A = VGPrice(1369.41, 1400,0.50215,-0.22898,0.20722,0.0541,0.012,0.56164,-1)
#print(A)



