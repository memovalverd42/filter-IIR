import numpy as n



############################## CONFIGURACION ######################################

Ap  = 0.5               # Rizo en la banda pasante
As  = 35                # Atenuacion
F   = 600               # Frecuencia de muestreo

# Para LOWPASS y HIGHPASS
fp = 60                 # Frecuencia de polo
fs = 120                # Frecuencia de cero

# Para BANDPASS y BANDSTOP
fp1 = 60                # Frecuencia en el polo 1 
fp2 = 240               # Frecuencia en el polo 2 
fs1 = 120               # Frecuencia en el cero 1 
fs2 = 180               # Frecuencia en el cero 2 

F_n = 0                 # Frecuencia para normalizar el filtro

# filter_type   ->  'LP' | 'HP'
filter_type = 'HP'
# filter        ->  'B'  |  'C'  |   'E' 
filter      = 'B'

BP = False              # Para BANDPASS
BS = True               # Para BANDSTOP

def truncar(num, n):                                                        # Función para truncar valores
    truncado = int(num * (10**n))/(10**n)
    return float(truncado)

if BP and filter_type == 'HP':
    fs = fs1
    fp = fp1
elif BP and filter_type == 'LP':
    fs = fs2
    fp = fp2
elif BS and filter_type == 'LP':
    fs = fs1
    fp = fp1
    F_n = fp1/2
elif BS and filter_type == 'HP':
    fs = fs2
    fp = fp2
    F_n = fp2 + (((F/2)-fp2)/2)
elif BP and filter_type == 'LP':
    pass
elif BP and filter_type == 'HP':
    pass

###############################################################################################################################################

########## 1. Determinar orden del filtro ##########

def determinar_orden(K):                                                    # Funcion  para determinar el orden del filtro
    
    # Determinar N segun el filtro:
    if filter == 'B':
        N = (n.log10(A))/((n.log10(1/K))) if BS == False else (n.log10(A))/((n.log10(K)))                                  # Orden para Butterworth (Formula)
        return int(N+1)                                         
    elif filter == 'C':                                                    
        N = (n.arccosh(A))/(n.arccosh(1/K))                                 # Orden para Chevysheb (Formula)
        return int(N+1)
    elif filter == 'H':
        q0 = (1-((1-K**2)**0.25))/(2*(1+(1-K**2)**0.25))
        q = q0+(2*q0**5)+(15*q0**9)+(150*q0*13)

        N = (n.log10(16*(A**2)))/(n.log10(1/q))                             # Orden para Eliptico (Formula)  
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
    print(f'\nK = {K} | N = {Nr} | k = {k}\n')

########## 2. Determinar los polos en el plano S mediante el polinomio de Butterworth ##########

def polos():                                                                # Funcion para calcular los polos en plano S
    S_real = []
    S_imag = []
    for i in range(1, k+1):
        valueR = round(-n.sin((((2*i)-1)/(2*Nr))*n.pi), 10)
        valueI = round(n.cos((((2*i)-1)/(2*Nr))*n.pi), 10)
        # valueI = 0 if valueI <= 1.0e-10 else valueI

        S_real.append(valueR)
        S_imag.append(valueI)
    
    return S_real, S_imag

if filter == 'B':
    sr, si = polos()

print('\n*************** Polos en Plano S ******************')
for i in range(k):
    print(f'S{i+1} = {sr[i]} + - j{si[i]}')

########## 3. Transformar los polos del plano S al plano Z ##########

# Determinar parametro alpha

if filter_type == 'LP' and filter == 'B':
    alpha = ((10**(0.1*Ap)-1)**(-(1)/(2*Nr)))*n.tan((fp*n.pi)/(F))                                     # Butterworth LOWPASS parametro alpha
elif filter_type == 'HP' and filter == 'B':
    alpha = ((10**(0.1*Ap)-1)**(1/(2*Nr)))*n.tan((fp*n.pi)/(F))                                        # Butterworth HIGHPASS parametro alpha
elif filter == 'C':
    alpha = n.tan((n.pi*fp)/F)                                                                         # Chebyshev parametro alpha
elif filter == 'E':
    alpha = (n.tan((n.pi*fp)/F)*n.tan((n.pi*fp)/F))**0.5                                               # Eliptico parametro alpha


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

# for i in range(k):
#     print(term_2aZ[i])

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
            H = (Z+1)/(Z+(term_2aZ[i]/2))
            print(f'|H{i+1}(Z)| = {abs(H)} -> C{i+1} = {1/abs(H)}')
        elif filter_type == 'HP':
            H = (Z-1)/(Z+(term_2aZ[i]/2))
            print(f'|H{i+1}(Z)| = {abs(H)} -> C{i+1} = {1/abs(H)}')
        C.append(1/abs(H))
    else:
        # if term_2aZ[i] > 0:
        #     den = (Z_cuadrada+(term_2aZ[i]*Z)+term_a2b2[i])
        # else:
        #     den = (Z_cuadrada+(term_2aZ[i]*Z)+term_a2b2[i])

        if filter_type == 'LP':
            H = (Z_cuadrada+(2*Z)+1)/(Z_cuadrada+(term_2aZ[i]*Z)+term_a2b2[i])
        elif filter_type == 'HP':
            H = (Z_cuadrada-(2*Z)+1)/(Z_cuadrada+(term_2aZ[i]*Z)+term_a2b2[i])

        print(f'|H{i+1}(Z)| = {abs(H)} -> C{i+1} = {1/abs(H)}')
        C.append(1/abs(H))
print('\n')

# PARA SCRIPT DE MATLAB

for i in range(k):                                                      # Para listar coeficientes de normalizacion
    print(f'C{i+1}={C[i]};')
print('\n')

for i in range(k):                                                      # Para listar los denominadores
    if Zi[i].imag == 0:
        print(f'A{i+1}=[1, {term_2aZ[i]/2}];')
    else:
        print(f'A{i+1}=[1, {term_2aZ[i]}, {term_a2b2[i]}];')

print('\n')







