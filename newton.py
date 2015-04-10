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
        curX = curX - (fEval(curX) / fprimeEval(curX))

    return curX


def main():
    # This makes the "variable" for sympy
    x = Symbol('x')

    # The function we want to analyze
    y = x**3 - x**2 - x + 1.1
    
    root = newtonsMethod(y, x, 20.00)

    print "Final answer: ", root


if __name__ == '__main__':
    main()
