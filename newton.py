from sympy import *
import numpy as np
import scipy.stats
import math

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
            history = []
        numer = 2*fEval(curX)*fprimeEval(curX)
        denom = 2*fprimeEval(curX)*fprimeEval(curX) - fEval(curX)*fdoublePrimeEval(curX) 
        if denom == 0:
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
            return x2
        x1 = x2
        x2 = x2 - (fEval(x2) / slope)

    history.append(abs(fEval(x2)))
    return history

def convergenceRate(path):
    np.polyfit(range(len(path)),  np.log10(path), 1)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(range(len(path)), np.log10(path))
    return [slope, r_value**2]

def main():
    # This makes the "variable" for sympy
    x = Symbol('x')

    # The function we want to analyze
    y = x**3 - x**2 - x + 1.1
    
    path = secantMethod(y, x, 2.00, 2.10, .00000001)
    print "Final answer: ", path
    print convergenceRate(path)

    #root = secantMethod(y, x, 2.00, 3.00)
    #print "Final answer: ", root

if __name__ == '__main__':
    main()
