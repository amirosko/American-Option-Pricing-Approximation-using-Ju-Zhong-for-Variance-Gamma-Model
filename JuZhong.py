from BlackScholesEuropean import EuropeanOption
from scipy.stats import norm
from math import log, sqrt, exp

def JuZhongPrice(S, K, r, q, sigma, T, phi):
    '''
    American option calculation based on Ju-Zhong paper
    '''
    if phi == 1 and q == 0:
        return EuropeanOption(S, K, r, q, sigma, T, phi, greekCal = True)
    alpha = 2 * r/(sigma**2)
    beta = 2 * (r-q)/(sigma**2)
    hTau = 1 - exp(-r*T)
    lambdaH = (-(beta-1) + phi * sqrt((beta-1)**2 + 4 * alpha/hTau))/2

    qInfty =  (1 - beta + phi * sqrt((beta - 1)**2 + 4*alpha))/2
    sInfty = K/(1 - 1/qInfty)
    hi = (-phi*(r-q)*T - 2*sigma*sqrt(T)) * K / (phi * (sInfty - K))
    initialGuess = sInfty + (K - sInfty) * exp(hi)
    Sx = findSx(initialGuess, K, r, q, sigma, T, phi, lambdaH)
    ah = (phi * (Sx - K) - EuropeanOption(Sx, K, r, q, sigma, T, phi=phi)[0])/hTau

    theta = EuropeanOption(Sx, K, r, q, sigma, T, phi=phi, greekCal=True)[-2]
    lambdaHDerivation = -phi * alpha / (hTau**2 * sqrt((beta -1)**2 + 4*alpha/hTau))
    b = (1 - hTau) * alpha * lambdaHDerivation/(2*(2 * lambdaH + beta - 1))
    c = - (1 - hTau) * alpha / (2 * lambdaH + beta - 1) * (-theta/(hTau * ah * r * exp(-r*T)) + 1/hTau + lambdaHDerivation/(2 * lambdaH + beta - 1))
    euroPrice = EuropeanOption(S, K, r, q, sigma, T, phi=phi)[0]
    if phi * (Sx - S) > 0:
        amerPrice = euroPrice + (hTau * ah * (S/Sx)**lambdaH)/(1 - b * (log(S/Sx))**2 - c * log(S/Sx))
    else:
        amerPrice = phi * (S - K)
    amerDelta, amerGamma, amerTheta, amerVega = greeksAnalysis(S, K, r, q, sigma, T, phi, Sx, ah, lambdaH, beta, alpha, hTau, b, c, amerPrice)
    amerPrice = float('%.3f' % amerPrice)
    amerDelta = float('%.3f' % amerDelta)
    amerGamma = float('%.3f' % amerGamma)
    amerVega = float('%.3f' % amerVega)
    amerTheta = float('%.3f' % amerTheta)
    euroPrice = float('%.3f' % euroPrice)
    return amerPrice, amerDelta, amerGamma, amerTheta, amerVega, amerPrice - euroPrice

def findSx(initialGuess, K, r, q, sigma, T, phi, lambdaH):
    finish = False
    countCycle = 0
    while (not finish):
        d1 = (log(initialGuess / K) + (r - q + sigma ** 2 / 2) * T) / (sigma * sqrt(T))
        euroPrice, delta = EuropeanOption(initialGuess, K, r, q, sigma, T, phi=phi, greekCal=True)[0:2]
        leftSide = phi * initialGuess - lambdaH * (phi * (initialGuess - K))
        rightSide = phi * initialGuess * exp(-q*T) * norm.cdf(phi*d1) - euroPrice * lambdaH
        # difference = phi * exp(-q*T) * norm.cdf(phi*d1) + lambdaH * (phi*(initialGuess-K) - euroPrice)/initialGuess - phi
        if abs(leftSide - rightSide) / K < 0.0000001:
            finish = True
        else:
            slopeBi = exp(-q*T) * norm.pdf(phi*d1) /(sigma*sqrt(T)) + (1 - lambdaH) * delta
            initialGuess = (lambdaH * K * phi + initialGuess * slopeBi - rightSide)/(slopeBi - phi * (1-lambdaH))
            countCycle += 1
    return initialGuess

def greeksAnalysis(S, K, r, q, sigma, T, phi, Sx, A, lambdaH, beta, alpha, h, b, c, amerPrice):
    euroPriceS, deltaS, gammaS, thetaS, vegaS = EuropeanOption(S, K, r, q, sigma, T, phi, greekCal=True)
    if phi * (Sx - S) > 0:
        euroPriceSx, deltaSx, gammaSx, thetaSx, vegaSx = EuropeanOption(Sx, K, r, q, sigma, T, phi, greekCal=True)
        d1Sx = (log(Sx / K) + (r - q + sigma ** 2 / 2) * T) / (sigma * sqrt(T))
        d2Sx = d1Sx - sigma * sqrt(T)
        d1S = (log(S / K) + (r - q + sigma ** 2 / 2) * T) / (sigma * sqrt(T))
        chi = b * (log(S/Sx))**2 + c * log(S/Sx)
        chiPS = (2 * b / S) * log(S/Sx) + c / S
        chiPSS = (2 * b / (S**2)) * (1 - log(S/Sx)) - c / (S**2)
        amerDelta = deltaS + (lambdaH/(S*(1-chi)) + chiPS/((1-chi)**2)) * (phi * (Sx - K) - euroPriceSx) * ((S/Sx)**lambdaH)
        amerGamma = exp(-q*T) * norm.pdf(phi*d1S) / (S * sigma * sqrt(T)) + (2 * lambdaH * chiPS / (S*(1-chi)**2) + 2 * chiPS**2 / ((1-chi)**3) + chiPSS / (1-chi)**2 + (lambdaH**2 - lambdaH)/(S**2 * (1-chi))) * (phi * (Sx - K) - euroPriceSx) * ((S/Sx)**lambdaH)
        amerTheta = r * amerPrice - (sigma*S)**2 * amerGamma / 2 - (r-q) * S * amerDelta
        # Vega
        '''
        paramLambda = (beta-1)**2 + 4 * alpha / h
        alphaP = - 2 * alpha / sigma
        betaP = - 2 * beta / sigma
        lambdaP = 0.5 * (-betaP + 0.5 * phi * (2*(beta-1)*betaP + 4 * alphaP/h)/(sqrt(paramLambda)))
        lambdaPh = - phi * alpha / (h**2 * sqrt(paramLambda))
        lambdaPhP = - (phi / (h**2)) * (alphaP/sqrt(paramLambda) - (alpha/2)*(2*(beta-1)*betaP + 4 * alphaP/h)/(paramLambda)**(3/2))
        bP = ((1-h)/2) * ((lambdaPh * alphaP + alpha * lambdaPhP) / (2*lambdaH + beta - 1) - alpha * lambdaPh * (2 * lambdaP + betaP) / (2*lambdaH + beta - 1)**2)
        numeratorForSxP = phi * (Sx - K) * lambdaP + Sx * exp(-q*T) * norm.pdf(d1Sx) * (sqrt(T) - d1Sx/sigma) - euroPriceSx * lambdaP - lambdaH * vegaSx
        denominatorSxP = phi - phi * lambdaH - deltaSx - exp(-q*T) * norm.pdf(d1Sx) / (sigma * sqrt(T)) + lambdaH * deltaSx
        SxP = numeratorForSxP / denominatorSxP
        d1SxP = SxP / (Sx * sigma * sqrt(T)) + (sqrt(T) - d1Sx / sigma)
        d2SxP = d1SxP - sqrt(T)
        AP = ((phi - deltaSx) * SxP - vegaSx) / h
        thetaP = - (exp(-q*T)/(2*sqrt(T))) * (Sx * norm.pdf(d1Sx) + sigma * norm.pdf(d1Sx) * SxP - Sx * sigma * norm.pdf(d1Sx) * d1Sx * d1SxP) - r * K * exp(-r*T) * norm.pdf(d2Sx) * d2SxP + q * Sx * exp(-q*T) * norm.pdf(d1Sx) * d1SxP + q * deltaSx * SxP
        cP = - (1-h) / (2 * lambdaH + beta - 1) * (alphaP - alpha * (2 * lambdaP + betaP) / (2 * lambdaH + beta - 1)) * (-thetaSx / (h * r * exp(-r*T) * A) + 1/h + lambdaPh / (2 * lambdaH + beta - 1)) - (1-h) * alpha / (2 * lambdaH + beta - 1) * ((-1 / (h * r * exp(-r*T))) * (thetaP / A - thetaSx * AP / (A**2)) + lambdaPhP / (2 * lambdaH + beta - 1) - lambdaPh * (2 * lambdaP + betaP) / ((2 * lambdaH + beta - 1)**2))
        amerVega = vegaS + (h/(1-chi)) * (AP + A * (lambdaP * log(S/Sx) - lambdaH * SxP / Sx)) * ((S/Sx)**lambdaH) - ((h*A*(S/Sx)**lambdaH)/(1-chi)**2) * (-bP * (log(S/Sx))**2 + 2 * b * log(S/Sx) * SxP / Sx - cP * log(S/Sx) + c * SxP / Sx)
        '''
        paramLambda = (beta - 1) ** 2 + 4 * alpha / h
        alphaP = - 2 * alpha / sigma
        betaP = - 2 * beta / sigma
        lambdaP = 0.5 * (-betaP + 0.5 * phi * (2 * (beta - 1) * betaP + 4 * alphaP / h) / sqrt(paramLambda))
        lambdaPh = - phi * alpha / (h**2 * sqrt(paramLambda))
        lambdaPhP = - (phi / ((h**2) * sqrt(paramLambda))) * (alphaP - alpha * (2 * (beta - 1) * betaP + 4 * alphaP / h) / (2 * paramLambda))
        bP = (1-h) * (lambdaPh * alphaP + alpha * lambdaPhP - alpha * lambdaPh * (2 * lambdaP + betaP) / (2 * lambdaH + beta - 1)) / (2 * (2 * lambdaH + beta - 1))
        numeratorForSxP = Sx * exp(-q*T) * norm.pdf(d1Sx) * (sqrt(T) - d1Sx / sigma) + (phi * (Sx - K) - euroPriceSx) * lambdaP - lambdaH * vegaSx
        denominatorSxP = (1 - lambdaH) * (phi - deltaSx) - exp(-q*T) * norm.pdf(d1Sx) / (sigma * sqrt(T))
        SxP = numeratorForSxP / denominatorSxP
        d1SxP = SxP / (Sx * sigma * sqrt(T)) + (sqrt(T) - d1Sx / sigma)
        d2SxP = d1SxP - sqrt(T)
        AP = (SxP * (phi - deltaSx) - vegaSx) / h
        thetaP = - exp(-q*T) * norm.pdf(d1Sx) * (Sx + sigma * SxP - Sx * sigma * d1Sx * d1SxP) / (2 * sqrt(T)) - r * K * exp(-r*T) * norm.pdf(d2Sx) * d2SxP + q * Sx * exp(-q*T) * norm.pdf(d1Sx) * d1SxP + q * deltaSx * SxP
        cP = - ((1-h) * alpha / (2 * lambdaH + beta - 1)) * (-(thetaP / A - thetaSx * AP / (A**2)) / (h * r * exp(-r*T)) + lambdaPhP / (2 * lambdaH + beta - 1) - lambdaPh * (2 * lambdaP + betaP) / ((2 * lambdaH + beta - 1)**2)) - (1-h) * (-thetaSx / (h * r * exp(-r*T) * A) + 1 / h + lambdaPh / (2 * lambdaH + beta - 1)) * (alphaP / (2 * lambdaH + beta - 1) - alpha * (2 * lambdaP + betaP) / ((2 * lambdaH + beta - 1)**2))
        chiP = bP * ((log(S/Sx))**2) - (2 * b * (log(S/Sx)) + c) * SxP / Sx + cP * (log(S/Sx))
        amerVega = vegaS + (h * ((S/Sx)**lambdaH) / (1 - chi)) * (AP + A * (lambdaP * (log(S/Sx)) - lambdaH * SxP / Sx)) + h * A * ((S/Sx)**lambdaH) *chiP / ((1 - chi)**2)

    else:
        amerDelta = phi
        amerGamma = 0
        amerTheta = 0
        amerVega = 0
    return amerDelta, amerGamma, amerTheta, amerVega

