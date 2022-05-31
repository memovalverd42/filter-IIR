import numpy as np
import filtros

filtro = filtros.Filter('B', 'LP', 500, 750, F = 2500, Ap = 0.6, As = 40)


print(filtro.parameter_K())



