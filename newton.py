from sympy import *
import numpy as np
import math

def newtonsMethod(f, x, startPoint):
    ''' Take a function f of a variable (symbol) x. Run Newton's Method
    starting at startPoint '''
    
    curX = startPoint
    
    # Sympy does this for us
    fprime = f.diff(x)

    # The lambdified version allow us to directly compute values
    fEval = lambdify(x, f, 'numpy')
    fprimeEval = lambdify(x, fprime, 'numpy')

    # Run until we get very close to a root
    while abs(fEval(curX) - 0) > 0.001:
        print curX, fEval(curX)
        slope = fprimeEval(curX)
        if fprimeEval(curX) == 0:
            return curX
        curX = curX - (fEval(curX) / fprimeEval(curX))

    return curX

def halleysMethod(f, x, startPoint):
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

    # Run until we get very close to a root
    while abs(fEval(curX) - 0) > 0.001:
        print curX, fEval(curX)
        numer = 2*fEval(curX)*fprimeEval(curX)
        denom = 2*fprimeEval(curX)*fprimeEval(curX) - fEval(curX)*fdoublePrimeEval(curX) 
        if denom == 0:
            return curX
        curX = curX - (numer / denom)

    return curX

def secantMethod(f, x, startPoint1, startPoint2):
    ''' Take a function f of a variable (symbol) x. Run Secant Method
    starting at startPoints '''
    
    x1 = startPoint1
    x2 = startPoint2
    
    # The lambdified version allow us to directly compute values
    fEval = lambdify(x, f, 'numpy')

    # Run until we get very close to a root
    while abs(fEval(x2) - 0) > 0.001:
        print x1, x2, fEval(x2)
        slope = (fEval(x2)-fEval(x1))/(x2-x1)
        if slope == 0:
            return x2
        x1 = x2
        x2 = x2 - (fEval(x2) / slope)

    return x2

def main():
    # This makes the "variable" for sympy
    x = Symbol('x')

    # The function we want to analyze
    y = x**3 - x**2 - x + 1.1
    
    root = halleysMethod(y, x, 2.00)
    print "Final answer: ", root

    #root = secantMethod(y, x, 2.00, 3.00)
    #print "Final answer: ", root

if __name__ == '__main__':
    main()
