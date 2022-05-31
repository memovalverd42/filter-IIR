import numpy as n

############################## CONFIGURACION ######################################

Ap  = 0.1737              # Rizo en la banda pasante
As  = 40               # Atenuacion
F   = 2500           # Frecuencia de muestreo

# Para LOWPASS y HIGHPASS
fp = 500                # Frecuencia de polo
fs = 750                # Frecuencia de cero

# Para BANDPASS y BANDSTOP
fp1 = 120                # Frecuencia en el polo 1 
fp2 = 180               # Frecuencia en el polo 2 
fs1 = 60                # Frecuencia en el cero 1 
fs2 = 240               # Frecuencia en el cero 2 

F_n = 0                    # Frecuencia para normalizar el filtro

# filter_type   ->  'LP' | 'HP'
filter_type = 'LP'
# filter        ->  'B'  |  'C'  |   'E' 
filter      = 'E'

BP = False              # Para BANDPASS
BS = False               # Para BANDSTOP

def truncar(num, n):                                                        # Función para truncar valores
    truncado = int(num * (10**n))/(10**n)
    return float(truncado)

if BP and filter_type == 'HP':
    fs = fs1
    fp = fp1
    F_n = fp1+((fp2-fp1)/2)
elif BP and filter_type == 'LP':
    fs = fs2
    fp = fp2
    F_n = fp1+((fp2-fp1)/2)
elif BS and filter_type == 'LP':
    fs = fs1
    fp = fp1
    F_n = fp1/2
elif BS and filter_type == 'HP':
    fs = fs2
    fp = fp2
    F_n = fp2 + (((F/2)-fp2)/2)
elif filter_type == 'LP':
    F_n = fp/2
elif filter_type == 'HP':
    F_n = fp + (((F/2)-fp)/2)

###############################################################################################################################################

########## 1. Determinar orden del filtro ##########
def determinar_orden(K):                                                    # Funcion  para determinar el orden del filtro
    
    # Determinar N segun el filtro:
    if filter == 'B':
        N = (n.log10(A))/((n.log10(1/K))) if BS == False else (n.log10(A))/((n.log10(K)))                                  # Orden para Butterworth (Formula)
        return int(N+1)                                         
    elif filter == 'C':                                                    
        N = (n.arccosh(A))/(n.arccosh(1/K)) if BS == False else (n.arccosh(A))/(n.arccosh(K))                               # Orden para Chevysheb (Formula)
        print(N)
        return int(N+1)
    elif filter == 'E':
        global q
        q0 = (1-((1-K**2)**0.25))/(2*(1+(1-K**2)**0.25))
        q = q0+(2*(q0**5))+(15*(q0**9))+(150*(q0**13))

        N = (n.log10(16*(A**2)))/(n.log10(1/q))                             # Orden para Eliptico (Formula)  
        print(f'\nq0 = {q0}     q = {q}')
        return int(N+1)

A = ((10**(0.1*As)-1)/(10**(0.1*Ap)-1))**0.5                                # A -> Formula estandar

# Determinar parametro K
if filter_type == 'LP' and BP == False and BS == False:
    K = (n.tan((n.pi*fp)/(F)))/(n.tan((n.pi*fs)/(F)))                     # K -> Formula Lowpass
    Nr = determinar_orden(K)
elif filter_type == 'HP' and BP == False and BS == False:
    K = (n.tan((n.pi*fs)/(F)))/(n.tan((n.pi*fp)/(F)))                     # K -> Formula Highpass
    Nr = determinar_orden(K)
if BP or BS:
    KA = n.tan((n.pi*fp2)/F)-n.tan((n.pi*fp1)/F)                            # Ka -> Formula 
    KB = n.tan((n.pi*fp1)/F)*n.tan((n.pi*fp2)/F)                            # Kb -> Formula 
    KC = n.tan((n.pi*fs1)/F)*n.tan((n.pi*fs2)/F)                            # Kc -> Formula 

    K1 = (KA*n.tan((n.pi*fs1)/F))/(KB-(n.tan((n.pi*fs1)/F))**2)             # K1 -> Formula 
    K2 = (KA*n.tan((n.pi*fs2)/F))/((n.tan((n.pi*fs2)/F))**(2)-KB)           # K2 -> Formula 

    if BP:                                                
        K = K1 if KC>=KB else K2                                            # K -> Para Bandpass
        Nr = determinar_orden(K)

    elif BS:                                   
        K = 1/K2 if KC>=KB else 1/K1                                        # K -> Para Bandstop
        Nr = determinar_orden(1/K)
    
k = int(Nr/2) if Nr%2 == 0 else int((Nr+1)/2)                               # Parametro k para el numero de iteraciones

if BP or BS:
    print(f'\nKA = {KA} | KB = {KB} | KC = {KC} | K1 = {K1} | K2 = {K2} | K = {K} | N = {Nr} | k = {k}\n | F_n = {F_n}')
else:
    print(f'\nA = {A} | K = {K} | N = {Nr} | k = {k}\n | F_n = {F_n}')

########## 2. Determinar los polos en el plano S mediante el polinomio de Butterworth ##########

# Determinar parametro alpha

if filter_type == 'LP' and filter == 'B':
    alpha = ((10**(0.1*Ap)-1)**(-(1)/(2*Nr)))*n.tan((fp*n.pi)/(F))                                     # Butterworth LOWPASS parametro alpha
elif filter_type == 'HP' and filter == 'B':
    alpha = ((10**(0.1*Ap)-1)**(1/(2*Nr)))*n.tan((fp*n.pi)/(F))                                        # Butterworth HIGHPASS parametro alpha
elif filter == 'C':
    alpha = n.tan((n.pi*fp)/F)                                                                         # Chebyshev parametro alpha
elif filter == 'E':
    alpha = (n.tan((n.pi*fp)/F)*n.tan((n.pi*fs)/F))**0.5                                               # Eliptico parametro alpha

def polos():                                                                # Funcion para calcular los polos en plano S
    S_real = []
    S_imag = []
    if filter == 'B':
        for i in range(1, k+1):
            valueR = round(-n.sin((((2*i)-1)/(2*Nr))*n.pi), 10)
            valueI = round(n.cos((((2*i)-1)/(2*Nr))*n.pi), 10)
            # valueI = 0 if valueI <= 1.0e-10 else valueI

            S_real.append(valueR)
            S_imag.append(valueI)

        return S_real, S_imag
    
    elif filter == 'C':
        epsi = (10**(0.1*Ap)-1)**0.5
        y    = (1/Nr)*n.arcsinh(1/epsi)
        x = [(2*i-1)*(n.pi/(2*Nr)) for i in range(1, k+1)]
        print(f'Epsilon = {epsi} | y = {y}\n')
        for i in range(1, k+1):
            valueR = round(-n.sin(x[i-1])*n.sinh(y), 10)
            valueI = round(n.cos(x[i-1])*n.cosh(y), 10)

            S_real.append(valueR)
            S_imag.append(valueI)

        return S_real, S_imag

    elif filter == 'E':
        lam = (1/(2*Nr))*(n.log(((10**(0.05*Ap))+1)/(((10**(0.05*Ap))-1))))
        sum_num = 0  
        sum_den = 0  
        for m in range(4):
            sum_num += ((-1)**m) * (q**(m*(m+1))) * (n.sinh(((2*m)+1)*lam))

        for m in range(1, 5):
            sum_den += ((-1)**m) * (q**(m**2)) * (n.cosh((2*m)*lam))

        ro = (-2*(q**0.25)*sum_num)/(1+2*sum_den)

        omega = []
        for i in range(1, k+1):
            u = i-0.5 if Nr%2 == 0 else i
            num = 0
            den = 0
            for m in range(4):
                num += ((-1)**m) * (q**(m*(m+1))) * (n.sin(((n.pi*u)/Nr)*((2*m)+1)))
            for m in range(1, 5):
                den += ((-1)**m) * (q**(m**2)) * (n.cos(((n.pi*u)/Nr)*(2*m)))

            omega.append((2*(q**0.25)*num)/(1+2*den))

        W = ( (1+(K*ro**2)) * (1+((ro**2)/K)) )**0.5

        V = [( (1-(K*omega[i]**2)) * (1-((omega[i]**2)/K)) )**0.5 for i in range(k)]

        for i in range(k):
            den_polos = 1+((ro**2)*(omega[i]**2))
            valueR = round((ro*V[i])/den_polos, 10)
            valueI = round((omega[i]*W)/den_polos, 10)

            S_real.append(valueR)
            S_imag.append(valueI)

        zeros = [(complex(0,1)/omega[i]).imag for i in range(k)]
        Zz_real = [( 1-( (zeros[i]**2)*(alpha**2) ) )/( 1+( (zeros[i]**2)*(alpha**2) ) ) for i in range(k)]
        Zz_img = [( 2*zeros[i]*alpha )/( 1+( (zeros[i]**2)*(alpha**2) ) ) for i in range(k)]

        # print(Zz_real)
        # print(Zz_img)
        # print(zeros)

        print(f'\nLambda = {lam} |  ro = {ro} | u = {u} | omega = {omega} | W = {W}')
        print(f'\n W = {W} | V = {V}')
        

        return S_real, S_imag, Zz_real, Zz_img, zeros

# if filter == 'B':
#     sr, si = polos()
# elif filter == 'C':
if filter != 'E': sr, si = polos() 
else: sr, si, Zz_r, Zz_i, zeros = polos()


print('\n*************** Polos en Plano S ******************')
for i in range(k):
    print(f'S{i+1} = {sr[i]} + - j{si[i]}')

if filter == 'E':
    print('\n*************** Zeros en Plano S ******************')
    for i in range(k):
        print(f'Sz{i+1} = j{zeros[i]}')

########## 3. Transformar los polos del plano S al plano Z ##########

def S_to_Z():
    if filter_type == 'LP':
        Di  = [1+(2*abs(sr[i])*alpha)+((sr[i])**2+(si[i])**2)*(alpha)**2 for i in range(k)]
        Zir = [(1-((sr[i])**2+(si[i])**2)*(alpha)**2)/(Di[i]) for i in range(k)]
        Zii = [(2*si[i]*alpha)/(Di[i]) for i in range(k)]
        Zi  = [complex(Zir[i], Zii[i]) for i in range(k)]

    elif filter_type == 'HP':
        Di  = [alpha**2+(2*abs(sr[i])*alpha)+((sr[i])**2)+((si[i])**2) for i in range(k)]
        Zir = [-((alpha**2)-((sr[i])**2+(si[i])**2))/(Di[i]) for i in range(k)]
        Zii = [(2*si[i]*alpha)/(Di[i]) for i in range(k)]
        Zi  = [complex(Zir[i], Zii[i]) for i in range(k)]

    term_2aZ = []
    term_a2b2 = []

    for i in range(k):
        if Zir[i] < 0:
            term_2aZ.append(2*abs(Zi[i].real))
        else:
            term_2aZ.append(-2*Zi[i].real)
        term_a2b2.append((Zi[i].real**2)+(Zi[i].imag**2))


    return Di, Zi, term_2aZ, term_a2b2

print(f'\nA = {A} | alpha = {alpha}\n')

Di, Zi, term_2aZ, term_a2b2 = S_to_Z()
if filter == 'E':
    term_2aZ_z  = [] 
    term_a2b2_z = []
    for i in range(k):
        if Zz_r[i] < 0:
            term_2aZ_z.append(2*abs(Zz_r[i]))
        else:
            term_2aZ_z.append(-2*Zz_r[i])
        term_a2b2_z.append((Zz_r[i]**2)+(Zz_i[i]**2))


print('\n*************** Di ******************')
for i in range(k):
    print(f'Di({i+1}) = {truncar(Di[i], 7)}')

print('\n*************** Zi ******************')
for i in range(k):
    print(f'Zi({i+1}) = {truncar(Zi[i].real, 7)} + - j{truncar(Zi[i].imag, 7)}')

print('\n*************** Denominador de H(Z) ******************')
for i in range(k):
    if Zi[i].imag == 0:
        print(f'H{i+1}(Z) = Z  {truncar(term_2aZ[i]/2, 7)}')
    else:
        print(f'H{i+1}(Z) = Z^2 {truncar(term_2aZ[i], 7)}Z + {truncar(term_a2b2[i], 7)}')


if filter == 'E':
    print('\n*************** Zeros en Plano Z ******************')
    for i in range(k):
        print(f'Zz{i+1} = {Zz_r[i]} + - j{Zz_i[i]}')

    print('\n*************** Numerador de H(Z) ******************')
    for i in range(k):
        if Zz_r[i] == 0:
            print(f'H{i+1}(Z) = Z  {truncar(term_2aZ[i]/2, 7)}') #PENDIENTE
        else:
            print(f'Hz{i+1}(Z) = Z^2 {truncar(term_2aZ_z[i], 7)}Z + {truncar(term_a2b2_z[i], 7)}')    

########## 4. Normalizar el filtro, tener ganancia unitaria ##########

Z = complex(round(n.cos((2*n.pi*F_n)/F), 8), n.sin((2*n.pi*F_n)/F))
Z_cuadrada = complex(n.cos((4*n.pi*F_n)/F), round(n.sin((4*n.pi*F_n)/F), 8))
print(f'\nZ = {Z}')
print(f'Z^2 = {Z_cuadrada}')

# print((Z_cuadrada+(2*Z)+1))

print('\n*************** H(Z) NORMALIZADA ******************')
C = []
for i in range(k):
    if Zi[i].imag == 0:
        if filter_type == 'LP':
            if filter != 'E':
                H = (Z+1)/(Z+(term_2aZ[i]/2))
                print(f'|H{i+1}(Z)| = {abs(H)} -> C{i+1} = {1/abs(H)}')
            else:
                H = (Z+1)/(Z+(term_2aZ[i]/2)) #PENDIENTE

        elif filter_type == 'HP':
            if filter != 'E':
                H = (Z-1)/(Z+(term_2aZ[i]/2))
                print(f'|H{i+1}(Z)| = {abs(H)} -> C{i+1} = {1/abs(H)}')
            else:
                H = (Z-1)/(Z+(term_2aZ[i]/2)) #PENDIENTE
        C.append(1/abs(H))
    else:
        if filter_type == 'LP':
            if filter != 'E': 
                H = (Z_cuadrada+(2*Z)+1)/(Z_cuadrada+(term_2aZ[i]*Z)+term_a2b2[i])
            else: 
                H = (Z_cuadrada+(term_2aZ_z[i]*Z)+term_a2b2_z)/(Z_cuadrada+(term_2aZ[i]*Z)+term_a2b2[i]) 
        elif filter_type == 'HP':
            if filter != 'E': 
                H = (Z_cuadrada-(2*Z)+1)/(Z_cuadrada+(term_2aZ[i]*Z)+term_a2b2[i])
            else: 
                H = (Z_cuadrada+(term_2aZ_z*Z)+term_a2b2_z)/(Z_cuadrada+(term_2aZ[i]*Z)+term_a2b2[i])
        
        # print(H, type(H))
        print(f'|H{i+1}(Z)| = {abs(H)} -> C{i+1} = {1/abs(H)}')
        C.append(1/abs(H))
print('\n')

# PARA SCRIPT DE MATLAB

for i in range(k):                                                      # Para listar coeficientes de normalizacion
    if filter != 'E': print(f'C{i+1}={C[i]};')
    else: print(f'C{i+1}={C[i][i]};')

print('\n')

for i in range(k):                                                      # Para listar los denominadores
    if Zi[i].imag == 0:
        print(f'A{i+1}=[1, {term_2aZ[i]/2}];')
    else:
        print(f'A{i+1}=[1, {term_2aZ[i]}, {term_a2b2[i]}];')

print('\n')







