from functools import reduce
import numpy as np

# Multiplicacion (Z+A+jB)(Z-A+jB)
def multi_complex(n: complex, k):
    term_2aZ = term_a2b2 = []
    for i in range(k):
        if n[i].real < 0: term_2aZ.append(2*abs(n[i].real))
        else: term_2aZ.append(-2*n[i].real)
        term_a2b2.append( (n[i].real**2) + (n[i].imag**2) )

    return term_2aZ, term_a2b2