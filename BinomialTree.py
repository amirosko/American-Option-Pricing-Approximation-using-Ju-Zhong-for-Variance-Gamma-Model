from math import log, sqrt, exp
import numpy as np

def mainBinomial(S, K, r, q, sigma, T, phi, stepSize, sigmaDiff):
    '''
    Main function for binomial tree
    '''
    def doubleTest(S, K, r, q, sigma, T, phi, stepSize):
        result1 = AmerBinomialTree(S, K, r, q, sigma, T, phi, stepSize)
        result2 = AmerBinomialTree(S, K, r, q, sigma, T, phi, stepSize + 1)
        return [(x+y)/2 for x, y in zip(result1, result2)]

    amerPrice, amerDelta, amerGamma, amerTheta = doubleTest(S, K, r, q, sigma, T, phi, stepSize)
    lowerPrice = doubleTest(S, K, r, q, sigma-sigmaDiff, T, phi, stepSize)[0]
    upperPrice = doubleTest(S, K, r, q, sigma+sigmaDiff, T, phi, stepSize)[0]
    amerVega = (upperPrice - lowerPrice)/(2*sigmaDiff)
    amerPrice = float('%.3f' % amerPrice)
    amerDelta = float('%.3f' % amerDelta)
    amerGamma = float('%.3f' % amerGamma)
    amerVega = float('%.3f' % amerVega)
    amerTheta = float('%.3f' % amerTheta)
    return amerPrice, amerDelta, amerGamma, amerTheta, amerVega

def AmerBinomialTree(S, K, r, q, sigma, T, phi, stepSize):
    deltaT = T/stepSize
    stockPriceTree = np.zeros((stepSize+1, stepSize+1))
    stockPriceTree[0, 0] = S
    optionPriceTree = np.zeros((stepSize+1, stepSize+1))
    u = exp(sigma * sqrt(deltaT))
    d = 1/u
    rfRate = exp(r*deltaT)
    p = (rfRate * exp(-q*deltaT) - d)/(u-d)
    # Form stock tree
    for index in range(1, stepSize+1):
        stockPriceTree[index, index] = stockPriceTree[index-1, index-1] * u
    for i in range(stepSize + 1):
        for j in range(i+1, stepSize + 1):
            stockPriceTree[i, j] = stockPriceTree[i, j-1] * d
    # Construct Option Tree:
    optionPriceTree[:, -1] = np.maximum(phi * (stockPriceTree[:, -1] - K), 0)
    for j in range(stepSize-1, -1, -1):
        for i in range(j+1):
            val = (p * optionPriceTree[i+1,j+1] + (1-p) * optionPriceTree[i, j+1])/rfRate
            if val < phi * (stockPriceTree[i, j] - K):
                val = phi * (stockPriceTree[i, j] - K)
            optionPriceTree[i, j] = val
    # Greeks calculation
    ''' delta = (Cu-Cd)/(Su-Sd)'''
    delta = (optionPriceTree[1, 1] - optionPriceTree[0, 1])/(stockPriceTree[1, 1] - stockPriceTree[0, 1])
    '''delta1 = (Cuu-Cud)/(Suu-S), delta2 = (Cud-Cdd)/(S-Sdd), deltaH = 0.5(Suu+S) - 0.5(S-Sdd), gamma = (delta1-delta2)/deltaH'''
    gamma = ((optionPriceTree[2, 2]-optionPriceTree[1, 2])/(stockPriceTree[2, 2] - stockPriceTree[1, 2]) - (optionPriceTree[1, 2]-optionPriceTree[0, 2])/(stockPriceTree[1, 2] - stockPriceTree[0, 2]))/(0.5*(stockPriceTree[2, 2] - stockPriceTree[0, 2] ))
    '''theta = (Cud - C)(2 deltaT)'''
    theta = (optionPriceTree[1, 2]-optionPriceTree[0, 0])/(2*deltaT)
    return optionPriceTree[0, 0], delta, gamma, theta
