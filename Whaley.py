from BlackScholesEuropean import EuropeanOption
from scipy.stats import norm
from math import log, sqrt, exp

def BAWPrice(S, K, r, q, sigma, T, phi):
    '''
    Barone-Adesi and Whaley quadratic approximation for American vanilla options
    Finds the American prices as the European prices + premium
    Premium is based on Sx, the critical stock price above or below
    which it becomes optimal to exercise the American options
    C(S,K)=h(T)g(s,h)
    h = k = 1 - e^{-T}; alpha = M = 2r/sigma^2; beta = N = 2(r-q)/sigma^2
    lambdaH represents q1 q2 in paper
    :param S: Spot price
    :param K: Strike price
    :param r: Risk free rate
    :param q: Dividend yield
    :param sigma: Volatility
    :param T: Maturity
    :param phi: phi = 1: 'Call' ; phi = -1:  'Put'
    :return: American option price
    '''
    if phi == 1 and q == 0:
        return EuropeanOption(S, K, r, q, sigma, T, phi, greekCal = True)
    euroPrice = EuropeanOption(S, K, r, q, sigma, T, phi, greekCal = False)[0]
    beta = 2 * (r-q)/(sigma**2) # N in paper
    alpha = 2 * r / (sigma**2) # M
    h = 1 - exp(-r*T) # K in paper
    lambdaH = (1-beta + phi * sqrt((beta-1)**2 + 4*alpha/h))/2
    lammbdaInf = (1-beta + phi * sqrt((beta-1)**2 + 4*alpha))/2
    sxInfty = K / (1 - 1/lammbdaInf)
    hI = -((r-q)*T + 2*phi * sigma * sqrt(T)) * K / (sxInfty-K)
    initialGuess = sxInfty + (K - sxInfty) * exp(hI)
    Sx = findSx(initialGuess, K, r, q, sigma, T, lambdaH, phi)
    d1 = (log(Sx / K) + (r - q + sigma ** 2 / 2)*T) / (sigma * sqrt(T))
    A = phi * (Sx / lambdaH) * (1 - exp(-q*T) * norm.cdf(phi*d1))
    if phi*(Sx-S) > 0:
        amerPrice = euroPrice + A * (S/Sx)**lambdaH
    else:
        amerPrice = phi * (S-K)
    amerDelta, amerGamma, amerTheta, amerVega = greeksAnalysis(S, K, r, q, sigma, T, phi, Sx, A, d1, lambdaH, beta, alpha, h, amerPrice)
    amerPrice = float('%.3f' % amerPrice)
    amerDelta = float('%.3f' % amerDelta)
    amerGamma = float('%.3f' % amerGamma)
    amerVega = float('%.3f' % amerVega)
    amerTheta = float('%.3f' % amerTheta)
    euroPrice = float('%.3f' % euroPrice)
    return amerPrice, amerDelta, amerGamma, amerTheta, amerVega, amerPrice - euroPrice

def findSx(initialGuess, K, r, q, sigma, T, lambdaH, phi):
    '''
    Use  method in paper for estimating Sx (page 309)
    '''
    finish = False
    countCycle = 0
    while (not finish):
        d1 = (log(initialGuess / K) + (r - q + sigma ** 2 / 2) * T)/(sigma * sqrt(T))
        euroPrice, delta = EuropeanOption(initialGuess, K, r, q, sigma, T, phi, True)[0:2]
        leftSide = (initialGuess - K) * phi
        rightSide = euroPrice + phi * initialGuess * (1 - exp(-q*T) * norm.cdf(phi*d1)) / lambdaH
        if abs(leftSide - rightSide) / K < 0.0000001:
            finish = True
        else:
            slopeBi = delta * (1 - 1/lambdaH) + phi * (1 - phi * exp(-q*T) * norm.pdf(phi*d1)/(sigma * sqrt(T))) / lambdaH
            initialGuess = (K + phi * rightSide - phi * slopeBi * initialGuess)/(1 - phi * slopeBi)
            countCycle += 1
    return initialGuess

def greeksAnalysis(S, K, r, q, sigma, T, phi, Sx, A, d1Sx, lambdaH, beta, alpha, h, amerPrice):
    euroPriceS, deltaS, gammaS, thetaS, vegaS = EuropeanOption(S, K, r, q, sigma, T, phi, greekCal=True)
    if phi * (Sx - S) > 0:
        euroPriceSx, deltaSx, gammaSx, thetaSx, vegaSx = EuropeanOption(Sx, K, r, q, sigma, T, phi, greekCal=True)
        amerDelta = deltaS + A * lambdaH * (S**(lambdaH-1)) / (Sx**lambdaH)
        amerGamma = gammaS + A * lambdaH * (lambdaH - 1) * (S**(lambdaH-2)) / (Sx**lambdaH)
        amerTheta = r * amerPrice - (sigma*S)**2 * amerGamma / 2 - (r-q) * S * amerDelta
        # Vega
        alphaP = - 2 * alpha / sigma
        betaP = - 2 * beta / sigma
        lambdaP = 0.5 * (-betaP + 0.5 * phi * (2*(beta-1)*betaP + 4 * alphaP/h)/(sqrt((beta-1)**2 + 4 * alpha / h)))
        numeratorForSxP = vegaSx - Sx * exp(-q*T) * norm.pdf(d1Sx) * (sqrt(T) - d1Sx/sigma) / lambdaH - (phi - deltaSx) * Sx * lambdaP / (lambdaH**2)
        denominatorSxP = phi - deltaSx + exp(-q*T) * norm.pdf(d1Sx) / (lambdaH*sigma*sqrt(T)) - (phi - deltaSx)/lambdaH
        SxP = numeratorForSxP / denominatorSxP
        d1SxP = SxP/(Sx * sigma * sqrt(T)) + sqrt(T) - d1Sx/sigma
        AP = (SxP / lambdaH - Sx * lambdaP / (lambdaH**2)) * (phi - deltaSx) - exp(-q*T) * norm.pdf(phi*d1Sx) * d1SxP * (Sx/lambdaH)
        amerVega = vegaS + ((S/Sx)**lambdaH) * (AP + A * (lambdaP * log(S/Sx) - lambdaH * SxP / Sx))
    else:
        amerDelta = phi
        amerGamma = 0
        amerTheta = 0
        amerVega = 0

    return amerDelta, amerGamma, amerTheta, amerVega

