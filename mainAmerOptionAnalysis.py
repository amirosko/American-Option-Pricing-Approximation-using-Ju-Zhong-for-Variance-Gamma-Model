from BinomialTree import mainBinomial
from JuZhong import JuZhongPrice
from Whaley import BAWPrice
from VG import VGPrice
import matplotlib.pyplot as plt

#def priceGreek(S, K, r, q, sigma, T, phi, stepSize, sigmaDiff):
#    # S, K, r, q, sigma, T, phi
#    a = JuZhongPrice(S, K, r, q, sigma, T, phi)[:5]
#    b = BAWPrice(S, K, r, q, sigma, T, phi)[:5]
#    # S, K, r, q, sigma, T, phi, stepSize, sigmaDiff
#    c = mainBinomial(S, K, r, q, sigma, T, phi, stepSize, sigmaDiff)
#    outPut = ''
#    for i in range(5):
#        outPut = outPut + str(c[i]) + ' & '
#        outPut =outPut + str(b[i]) + ' & '
#        outPut =outPut + str(a[i]) + ' & '
#    return(outPut)
#outPut = priceGreek(120,100, 0.07, 0.03, 0.3, 3, 1,1000,0.001)
#print(outPut)
#print('JuZhong:',a)
#print('Whaley:',b)
#print('True', c)


def earlyExecerisePremium(S0,r,q,sigma,nu,theta,T,phi, strikeList):
    outPutEuro = []
    outPutAmer = []
    outPutEarly = []
    for K in strikeList:
        earlyExecerise, euroValue = VGPrice(S0,K,nu,theta,sigma,r,q,T,phi)
        earlyExecerise = float('%.4f' % earlyExecerise)
        euroValue = float('%.4f' % euroValue)
        outPutEuro.append(euroValue)
        outPutAmer.append(earlyExecerise + euroValue)
        outPutEarly.append(earlyExecerise)
    return outPutEuro, outPutAmer, outPutEarly

    #percentage = [x/y for x, y in zip(outPut, trueValueList)]
    #print(outPut)
    #plt.plot(strikeList, outPut, 'o', label='Approximation')
    #plt.plot(strikeList, trueValueList, 'or', label='True value')
    #plt.plot(strikeList, percentage, 'o', label = 'Percentage capture')
    #plt.title('T = {}'.format(T))
    #plt.xlabel('K')
    #plt.ylabel('Early Exercise Premium')
    #plt.ylabel('Percentage capture')
    #plt.legend(loc='upper left')
    #plt.show()

K = [K for K in range(1200, 1420, 20)]
#trueVal4 = [1.041,1.25,1.489,1.742,2.023,2.339,2.741,3.178,3.687,4.241,4.861]
A4 = earlyExecerisePremium(1369.41,0.0541,0.012,0.20722,0.50215,-0.22898,0.56164,-1, K)[-1]
#trueVal3 = [0.539, 0.677, 0.835, 0.997, 1.185, 1.388, 1.681, 2.008, 2.392, 2.795, 3.291]
#A3 = earlyExecerisePremium(1369.41,0.0549,0.011,0.19071,0.49083,-0.28113,0.46575,-1, K)[-1]
#trueVal2 = [0.063, 0.124, 0.191, 0.253, 0.322, 0.392, 0.602, 0.791, 1.078, 1.374, 1.678]
#A2 = earlyExecerisePremium(1369.41,0.0536,0.012,0.185,0.2246,-0.28837,0.21643,-1, K)[-1]
#trueVal1 = [0.025, 0.070, 0.117,0.159,0.196,0.239, 0.428, 0.608, 0.869, 1.103, 1.483]
#A1 = earlyExecerisePremium(1369.41,0.0533,0.011,0.17875,0.13317,-0.30649,0.13972,-1, K)[-1]
print(A4)

