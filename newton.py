#!/usr/bin/env python

from sympy import *
import numpy as np
import scipy.stats
import math
import gc
import time

def recursiveMethod(f, x, topLevel, degree, startPoint, epsilon):
    print "Degree:", degree, "f:", f
    if degree < 2:
        ret = solve(f)[0].evalf()
        #print "Returning", ret
        return ret
        #return 
    
    fEval = lambdify(x, f, 'numpy')

    # Recurse on a polynomial of odd degree
    newDeg = degree / 2
    if newDeg % 2 == 0:
        newDeg -= 1

    # fApprox is the taylor expansion around the point z
    z = Symbol('z')
    taylor = series(f, x, z, newDeg + 1).removeO()
    taylorAt = lambdify(z, taylor, 'numpy')
    
    curGuess = startPoint
    history = []
    
    while abs(fEval(curGuess)) > epsilon:
        if topLevel:
            history.append(abs(fEval(curGuess)))
            if len(history) > 1 and history[-1] > history[-2]:
                history = []
        
        # Compute Taylor expansion for degree / 2
        curGuess = recursiveMethod(taylorAt(curGuess), x, False, newDeg, curGuess, epsilon)
    
            
    if topLevel:
        history.append(abs(fEval(curGuess)))
        print "Found root:", curGuess
        return history
    else:
        return curGuess
    

def newtonsMethod(f, x, startPoint, epsilon):
    ''' Take a function f of a variable (symbol) x. Run Newton's Method
    starting at startPoint, and return the final monotonically decreasing
    sequence of decreasing value evaluations when it was evaluated'''

    curX = startPoint

    # Sympy does this for us
    fprime = f.diff(x)

    # The lambdified version allow us to directly compute values
    fEval = lambdify(x, f, 'numpy')
    fprimeEval = lambdify(x, fprime, 'numpy')

    history = []

    # Run until we get very close to a root
    while abs(fEval(curX) - 0) > epsilon:
        history.append(abs(fEval(curX)))
        if len(history) > 1 and history[-1] > history[-2]:
            history = []
        slope = fprimeEval(curX)
        if fprimeEval(curX) == 0:
            return []
        curX = curX - (fEval(curX) / fprimeEval(curX))

    history.append(abs(fEval(curX)))
    return history

def halleysMethod(f, x, startPoint, epsilon):
    ''' Take a function f of a variable (symbol) x. Run Halley's Method
    starting at startPoint '''

    curX = startPoint

    # Sympy does this for us
    fprime = f.diff(x)

    # Sympy does this for us
    fdoublePrime = fprime.diff(x)

    # The lambdified version allow us to directly compute values
    fEval = lambdify(x, f, 'numpy')
    fprimeEval = lambdify(x, fprime, 'numpy')
    fdoublePrimeEval = lambdify(x, fdoublePrime, 'numpy')
    history = []

    # Run until we get very close to a root
    while abs(fEval(curX) - 0) > epsilon:
        history.append(abs(fEval(curX)))
        if len(history) > 1 and history[-1] > history[-2]:
            history = [history[-1]]
        numer = 2*fEval(curX)*fprimeEval(curX)
        denom = 2*fprimeEval(curX)*fprimeEval(curX) - fEval(curX)*fdoublePrimeEval(curX)
        if denom == 0 or numer == 0:
            return []
        curX = curX - (numer / denom)

    history.append(abs(fEval(curX)))
    return history

def secantMethod(f, x, startPoint1, startPoint2, epsilon):
    ''' Take a function f of a variable (symbol) x. Run Secant Method
    starting at startPoints '''

    x1 = startPoint1
    x2 = startPoint2

    # The lambdified version allow us to directly compute values
    fEval = lambdify(x, f, 'numpy')

    history = []

    # Run until we get very close to a root
    while abs(fEval(x2) - 0) > epsilon:
        history.append(abs(fEval(x2)))
        if len(history) > 1 and history[-1] > history[-2]:
            history = []
        slope = (fEval(x2)-fEval(x1))/(x2-x1)
        if slope == 0:
            return []
        x1 = x2
        x2 = x2 - (fEval(x2) / slope)

    history.append(abs(fEval(x2)))
    return history

def convergenceRate(path):
    while path[-1] == 0.0:
        path = path[0:-1]
    logFArray = np.log10(path)
    #print logFArray
    #print logFArray[1:], " ", logFArray[0:-1]
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(logFArray[0:-1], logFArray[1:])
    return [slope, r_value**2]

def main():
    # This makes the "variable" for sympy
    x = Symbol('x')

    # The function we want to analyze
    y = x**9 - x**8 + x**7 - x**6 + x**5 - x**4 + x**3 - x**2 + x - 1

    #print "Recrusive method on (" + repr(y) + "):"
    #print recursiveMethod(y, x, True, 9, 5.0, 0.001)
    
        
    total = [0, 0]
    count = 0
    print "Recursive Method:"
    gc.disable()
    tic = time.time()
    for i in range(-1000, 1000, 100):
        path = recursiveMethod(y, x, True, 9, i/100., .01)
        path = [float(k) for k in path]
        if len(path) < 3:
            continue
        count = count+ 1
        total[0] = total[0] + convergenceRate(path)[0]
        total[1] = total[1] + convergenceRate(path)[1]
    toc = time.time()
    gc.enable()
    print total[0] / count, " ", total[1]/count
    print "AVG TIME: {}".format((toc-tic) / len(xrange(-1000, 1000, 100)))

#   
#     total = [0, 0]
#     count = 0
#     print "Secant Method:"
#     for i in range(-1000, 1000, 10):
#         path = secantMethod(y, x, i/100., (i+10)/100., .000000000001)
#         if len(path) < 3:
#             continue
#         count = count+ 1
#         total[0] = total[0] + convergenceRate(path)[0]
#         total[1] = total[1] + convergenceRate(path)[1]
#     print total[0] / count, " ", total[1]/count
# 
    total = [0, 0]
    count = 0
    print "Newton's Method:"
    gc.disable()
    tic = time.time()
    for i in range(-1000, 1000, 100):
        path = newtonsMethod(y, x, i/100., .01)
        if len(path) < 3:
            continue
        count = count+ 1
        total[0] = total[0] + convergenceRate(path)[0]
        total[1] = total[1] + convergenceRate(path)[1]
    toc = time.time()
    gc.enable()
    print total[0] / count, " ", total[1]/count
    print "AVG TIME: {}".format((toc-tic) / len(xrange(-1000, 1000, 100)))
# 
#     total = [0, 0]
#     count = 0
#     print "Halley's Method:"
#     for i in range(-1000, 1000, 10):
#         path = halleysMethod(y, x, i/100.,  .000000000001)
#         if len(path) < 3:
#             continue
#         count = count+ 1
#         total[0] = total[0] + convergenceRate(path)[0]
#         total[1] = total[1] + convergenceRate(path)[1]
#     print total[0] / count, " ", total[1]/count

if __name__ == '__main__':
    main()
