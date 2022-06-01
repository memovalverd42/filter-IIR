import numpy as np
from tools import operations

PI = np.pi

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
    def orden(self, K, A):
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
            return int(N+1)

    # Metodo para determinar parametro K del filtro
    def parameter_K_and_N(self):
        global K, N, k
        A = ( ( 10**(0.1*self.As)-1 )/( 10**(0.1*self.Ap)-1 ) )**0.5      # Parametro A
        if self.type == 'LP':
            K = ( np.tan( (PI*self.fp)/(self.F) ) )/( np.tan( (PI*self.fs)/(self.F) ) )  # Formula para K, Paso Bajo
            N = self.orden(K, A)
            # return K
        if self.type == 'HP':
            K = ( np.tan( (PI*self.fs)/(self.F) ) )/( np.tan( (PI*self.fp)/(self.F) ) )                     # K -> Formula Highpass
            N = self.orden(K, A)

        k = int( N/2 ) if N%2 == 0 else int( (N+1)/2 )

        return K, A, N, k 
    
    def constant_alpha(self):
        global alpha
        if self.type == 'LP' and self.filt == 'B':
            alpha = ( ( 10**(0.1*self.Ap)-1 )**( -(1)/(2*N) ) ) * np.tan( (self.fp*PI)/(self.F) )                                     # Butterworth LOWPASS parametro alpha
        elif self.type == 'HP' and self.filt == 'B':
            alpha = ( ( 10**(0.1*self.Ap)-1 )**( 1/(2*N) ) ) * np.tan( (self.fp*PI)/(self.F) )                                        # Butterworth HIGHPASS parametro alpha
        elif self.filt == 'C':
            alpha = np.tan( (PI*self.fp)/self.F )                                                                         # Chebyshev parametro alpha
        elif self.filt == 'E':
            alpha = ( np.tan( (PI*self.fp)/self.F ) * np.tan( (PI*self.fs)/self.F ) )**0.5 
        
        return alpha
    
    # Metodo para calcular polos en plano S
    def polos_S(self):
        global Ps
        if self.filt == 'B':
            Ps = [complex(round(-np.sin( ( ((2*i)-1 )/(2*N) ) * PI), 10), 
                        round(np.cos( ( ((2*i)-1)/(2*N) ) * PI ), 10)) 
                        for i in range(1, k+1)]
            
            return Ps

        elif self.filt == 'C':
            epsi = ( 10**( 0.1*self.Ap ) - 1 )**0.5
            y    = (1/N) * np.arcsinh( 1/epsi )
            x    = [(2*i-1) * ( PI/(2*N) ) for i in range(1, k+1)]

            Ps = [complex(round(-np.sin( x[i-1] ) * np.sinh(y), 10), 
                        round(np.cos( x[i-1] ) * np.cosh(y), 10)) 
                        for i in range(1, k+1)]
            
            return Ps, epsi, y, x

        elif self.filt == 'E':
            lam        = (1/(2*N)) * ( np.log( ((10**(0.05*self.Ap))+1) / (((10**(0.05*self.Ap))-1)) ) )
            sum_num_ro = np.sum([ ( (-1)**m ) * ( q**(m*(m+1)) ) * ( np.sinh( ((2*m)+1)*lam ) ) for m in range(4)])
            sum_den_ro = np.sum([ ( (-1)**m ) * ( q**(m**2) ) * ( np.cosh( (2*m)*lam ) ) for m in range(1, 5)])
            ro         = ( -2 * (q**0.25) * sum_num_ro ) / ( 1 + 2 * sum_den_ro )

            W = ( (1+(K*ro**2)) * (1+((ro**2)/K)) )**0.5

            global omega
            omega = []
            for i in range(1, k+1):
                u = i-0.5 if N%2 == 0 else i
                sum_num_omega = np.sum([ ((-1)**m) * (q**(m*(m+1))) * (np.sin(((PI*u)/N)*((2*m)+1))) for m in range(4)])
                sum_den_omega = np.sum([ ((-1)**m) * (q**(m**2)) * (np.cos(((PI*u)/N)*(2*m))) for m in range(1, 5)])

                omega.append( ( 2 * (q**0.25) * sum_num_omega ) / ( 1 + 2 * sum_den_omega ) )
            
            V = [( (1-(K*omega[i]**2)) * (1-((omega[i]**2)/K)) )**0.5 for i in range(k)]

            Ps = [complex(round( (ro*V[i]) / ( 1 + ( (ro**2) * (omega[i]**2) ) ), 10 ), 
                        round( (omega[i]*W) / ( 1 + ( (ro**2) * (omega[i]**2) ) ), 10 )) 
                        for i in range(k)]

            return Ps, lam, ro, W, omega, V

    # Metodo para calcular Zeros en el plano S (Caso de filtro Eliptico)
    def zeros_S(self):
        global Zs
        Zs = [ ( complex(0,1) / omega[i] ).imag for i in range(k) ]
        return Zs

    def S_to_Z(self):
        global Pz, Di
        if self.type == 'LP':
            Di  = [ 1 + (2*abs(Ps[i].real)*alpha) + ( (Ps[i].real)**2+(Ps[i].imag)**2 ) * (alpha)**2 for i in range(k) ]
            Pz = [complex(( 1 - ( (Ps[i].real)**2 + (Ps[i].imag)**2 ) * (alpha)**2 ) / (Di[i]), 
                            ( 2*Ps[i].imag*alpha ) / (Di[i])) 
                            for i in range(k)]

        elif self.type == 'HP':
            Di  = [alpha**2+(2*abs(Ps[i].real)*alpha)+((Ps[i].real)**2)+((Ps[i].imag)**2) for i in range(k)]
            Pz  = [ complex( -( (alpha**2) - ( (Ps[i].real)**2 + (Ps[i].imag)**2 ) ) / (Di[i]), 
                            ( 2*Ps[i].imag*alpha ) / (Di[i]) ) 
                            for i in range(k) ]

        if self.filt == 'E':
            Zz_real = [( 1-( (Zs[i]**2)*(alpha**2) ) )/( 1+( (Zs[i]**2)*(alpha**2) ) ) for i in range(k)]
            Zz_img = [( 2*Zs[i]*alpha )/( 1+( (Zs[i]**2)*(alpha**2) ) ) for i in range(k)]


        # Multiplicacion (Z+A+jB)(Z-A+jB)
        term_2aZ_p, term_a2b2_p = operations.multi_complex(Pz, k)
        # for i in range(k):
        #     if Pz[i].real < 0: term_2aZ.append(2*abs(Pz[i].real))
        #     else: term_2aZ.append(-2*Pz[i].real)
        #     term_a2b2.append( (Pz[i].real**2) + (Pz[i].imag**2) )

        if self.filt == 'E': 
            term_2aZ_z, term_a2b2_z = operations.multi_complex(Zz_real, k) 
            return Pz, Di, term_2aZ_p, term_a2b2_p, Zz_real, term_2aZ_z, term_a2b2_z


        return Pz, Di, term_2aZ_p, term_a2b2_p

    def normalization(self):
        pass
