# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 13:26:50 2023
@author: leasn
Python script to calculate the standard deviation. Does not return bits. 
"""

import math as m

Qi = [1152921504606584833, 35184372744193, 35184373006337, 35184368025601, 35184376545281] 
Pi = [2305843009211596801, 2305843009210023937]

logN = 16
h = 32768

depth = 4

levelQ = depth
levelP = 1

logScale = 45
logSlots = 15

std0 = 3.2
logN2 = logSlots/2 - 1
logN1 = 2 + logN2



    
#Write everything in a dictionary for better readability
params = {
    "logN":logN,
    "h": h,
    "depth" : depth,
    "logScale" : logScale,
    "logSlots" : logSlots,
    "Qi" : Qi ,
    "Pi" : Pi,
    "logN1" : logN1,
    "logN2" : logN2}



# productArray returns the product of the elements in the array
def productArray(array):
    q = 1
    for i in array:
        q *= i
    return q

# Number of sum during the gadget-product
def DecompRNS(levelQ, levelP):
    return (levelQ + levelP + 1) // (levelP + 1)

#Calculates the standard deviation of a fresh ciphertext
def stdFresh(params):
    return m.sqrt(std0**2*(2*params["h"] + 1))

#Calculates the standard deviation of the additional rounding error
def stdRound(params):
    return m.sqrt(1/12 + params["h"]/12.0)

def stdRescaling(params):
    pass

#Calculates the standard deviation of key switching
def stdKeySwitching(params, stdOld, levelQ):

    # Initial noise
    initial = stdOld**2

    # Inner Product Noise
    # (N * sigma**2 * sum(q_alpha_i)**2) / (12 * P**2) 

    sum_q_alphai = 0

    decompRNS = DecompRNS(levelQ, levelP)

    for i in range(decompRNS):

        start = i * (levelP+1)
        end = (i+1) * (levelP+1)

        if i == decompRNS-1:
            end = levelQ+1

        sum_q_alphai += productArray(Qi[start:end])

    inner_product = (1<<params["logN"]) * (sum_q_alphai**2) * (std0**2) / (12 * productArray(Pi)**2)

    # Rounding Noise
    rounding = stdRound(params)**2

    std = m.sqrt(initial + inner_product + rounding)
        
    return std

#Extract the real and the imaginary part. Only one return value, since the std are equal.
def stdExtractRealAndImag(params, stdOld, crtLevel):
    return stdOld**2 + 1/4*(stdKeySwitching(params, stdOld, crtLevel) - stdOld**2)

#Assembles the above for the std of CtS
def stdCoefficientsToSlots(params, stdOld, level):
    
    #Split the calculations in three factors, to make searching for errors easier.

    Q = productArray(Qi[:level+1]) #Are you sure ? Shouldn't it be Qi[level]? 
    
    #Compute first part of the sum
    tmpStD = stdOld**2 + stdKeySwitching(params, stdOld, level)**2
    
    boundCoeffs = (2**params["logScale"])**params["depth"]/(Q**2*(2**params["logSlots"]))
    
    rest = ((2**params["logN"])*(2**params["logSlots"]))**params["depth"]
    
    tmp = tmpStD*boundCoeffs*rest
    
    #Compute the ciphertext modulus
    for j in range(1, level+1):

        Q = productArray(Qi[:level+1-j]) #Are you sure ? Shouldn't it be Qi[level]? 

        #Compute the remaining sum
        boundCoeffs = (2**params["logScale"])**(params["depth"] - j)/(Q**2*(2**params["logSlots"]))
        
        rest = ((2**params["logN"])*(2**params["logSlots"]))**(params["depth"] - j)
        
        #Assemble all the summands.
        tmp += (stdKeySwitching(params, stdOld, level - j)**2)*boundCoeffs*rest        
    
    #Add final rounding error
    tmp += stdRound(params)**2
    
    return m.sqrt(stdExtractRealAndImag(params, m.sqrt(tmp), level))

#Assembles the above for the std of StC.
#Does not yet include the tweaking and the encoding of the matrices!
def stdSlotsToCoefficients(params, stdOld, crtLevel):
    tmp = 0
    for i in range(1,params["k"] + 1):
        tmp += (2**params["logN"]*2**params["logN1"]*2**params["logN2"])**i
    
    return m.sqrt(1.0/Qi[level]**2*(stdKeySwitching(params, stdOld, crtLevel)**2 - stdOld**2 - stdRound(params)**2)*tmp + stdRound(params)**2) 

def stdMatrixMul(params, stdOld, stdDiags, level):

    # |M|
    res = stdDiags**2

    # N * n1 * (old + ks) * |M|
    res *= (2**(params["logN"] + params["logN1"])) * (stdKeySwitching(params, stdOld, level)**2)
    
    # (N * n1 * (old + ks) * |M|) + 2round + ks - 2/12
    res += 2 * stdRound(params)**2
    res += stdKeySwitching(params, 0, level)**2
    res -= 2/12

    # n2 * (N * n1 * (old + ks) * |M|) + 2round + ks - 2/12
    res *= (2**params["logN2"])

    # (n2 * (N * n1 * (old + ks) * |M|) + 2round + ks - 2/12) + round
    res += stdRound(params)**2

    # (1/(q[level]**2)) * ((n2 * (N * n1 * (old + ks) * |M|) + 2round + ks - 2/12) + round)
    res /= Qi[level]**2

    # (1/(q[level]**2)) * ((n2 * (N * n1 * (old + ks) * |M|) + 2round + ks - 2/12) + round) + round
    res += stdRound(params)**2

    return m.sqrt(res)

def stdDFT(params, stdOld, stdDiags, level):

    if len(stdDiags) > level+1:
        raise Error("not enough level")

    for i in range(len(stdDiags)):
        stdOld = stdMatrixMul(params, stdOld, stdDiags[i], level-i)

    return stdOld

def stdC2S(params, stdOld, stdDiags, level, sparse:bool):
    
    #DFT
    dft = stdDFT(params, stdOld, stdDiags, level)

    #Conjugate
    dftConj = stdKeySwitching(params, dft, 0)

    #Real part extraction
    stdReal = m.sqrt((dft**2) + (dftConj**2))

    #If sparse: merges both real and imaginary
    if sparse:

        #Rotate imag by slots/2 (same std as real)
        stdImag = stdKeySwitching(params, stdReal, 0) 

        # Adds both ciphertexts
        stdReal = m.sqrt(stdReal**2 + stdImag**2)

    return stdReal


if __name__ == "__main__":
    
    #print(stdCoefficientsToSlots(params, stdFresh(params), levelQ))
    print(stdC2S(params, stdFresh(params), [1]*depth, levelQ, sparse=False))
