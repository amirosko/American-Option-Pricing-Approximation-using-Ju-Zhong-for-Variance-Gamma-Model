from scipy.stats import norm
from math import log, sqrt, exp
import matplotlib.pyplot as plt

def EuropeanOption(S, K, r, q, sigma, T, phi, greekCal = False):
    '''
    Calculation of Euro option price and greeks
    '''
    delta = None
    gamma = None
    theta = None
    vega = None

    top = log(S/K) + (r - q + sigma**2/2)*T
    bottom = sigma * sqrt(T)
    d1 = top/bottom
    d2 = d1 - sigma * sqrt(T)

    b1 = exp(-q*T)
    b2 = exp(-r*T)

    if greekCal:
        gamma = b1 * norm.pdf(d1)/(S * bottom)
        vega = b1 * S * norm.pdf(d1) * sqrt(T)

    if  phi == 1:
        nd1 = norm.cdf(d1)
        nd2 = norm.cdf(d2)
        price = S * b1 * nd1 - K * b2 * nd2
        if greekCal:
            delta = b1 * nd1
            theta = -b1 * S * norm.pdf(d1) * sigma / (2*sqrt(T)) - r * K * b2 * nd2 + q * S * b1 * nd1

    elif phi == -1:
        nNd1 = norm.cdf(-d1)
        nNd2 = norm.cdf(-d2)
        price = K * b2 * nNd2 - S * b1 * nNd1
        if greekCal:
            delta = -b1 * nNd1
            theta = -b1 * S * norm.pdf(d1) * sigma / (2*sqrt(T)) + r * K * b2 * nNd2 - q * S * b1 * nNd1

    return price, delta, gamma, theta, vega
