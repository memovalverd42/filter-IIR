import numpy as np
from tools import operations

class Filter:

    def __init__(self, filt: str, type: str, fp, fs, fp2 = 0, fs2 = 0, F = 0, Ap = 0.5, As = 40):
        self.filt   = filt
        self.type   = type
        self.fp     = fp
        self.fs     = fs
        self.fp2    = fp2
        self.fs2    = fp2
        self.F      = F
        self.Ap     = Ap
        self.As     = As
    
    # Metodo para determinar el orden del filtro
    def orden(self, K):
        A = ( ( 10**(0.1*self.As)-1 )/( 10**(0.1*self.Ap)-1 ) )**0.5      # Parametro A
        if self.filt == 'B':
            N = (np.log10(A))/((np.log10(1/K)))         # Orden para Butterworth (Formula)
            return int(N+1)
        elif self.filt == 'C':
            N = ( np.arccosh(A) )/( np.arccosh(1/K) )       # Orden para Chevysheb (Formula)
            return int(N+1)
        elif self.filt == 'E':
            global q
            q0 = ( 1-( (1-K**2)**0.25 ) )/( 2*( 1+(1-K**2)**0.25 ) )
            q  = q0 + ( 2*(q0**5) ) + ( 15*(q0**9) ) + ( 150*(q0**13) )

            N  = (np.log10(16*(A**2)))/(np.log10(1/q))                             # Orden para Eliptico (Formula)  
            print(f'\nq0 = {q0}     q = {q}') ## Empezar a usar format
            return int(N+1)

    # Metodo para determinar parametro K del filtro
    def parameter_K(self):
        global K, N
        # global N
        if self.type == 'LP':
            K = ( np.tan( (np.pi*self.fp)/(self.F) ) )/( np.tan( (np.pi*self.fs)/(self.F) ) )  # Formula para K, Paso Bajo
            N = self.orden(K)
            print(N)

    