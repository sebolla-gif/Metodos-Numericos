import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

X = np.array([-5,-3,-2,-1,0,2,4])
Y = np.array([-15,11,23,12,5,16,-5])

coeficientes_polinomio = np.polyfit(X, Y, len(X)-1)
x_fino = np.linspace(min(X), max(X), 1000)
Y_polinomio = np.polyval(coeficientes_polinomio, x_fino)
y1 = np.polyval(coeficientes_polinomio, 1)
print(f"El valor estimado de la funcion en 1 es: {y1}")

Spline = CubicSpline(X, Y)
x_fino = np.linspace(min(X), max(X), 1000)
y_spline = Spline(x_fino)
y1 = Spline(1)
print(f"El valor estimado de la funcion en 1 con spline es: {y1}")

plt.figure(figsize=(10, 4))
plt.axhline(0, color='k', linewidth=1)
plt.plot(X, Y, 'r*', markersize=8, label='Datos')
plt.plot(x_fino, Y_polinomio, 'b-', label='Polinomio interpolador')
plt.plot(x_fino, y_spline, 'g-', label='Polinomio interpolador Spline')
plt.title("Interpolación polinómica de los datos")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid()
plt.legend()
plt.show()

"""
En la tabla podemos observar que despues del x=0 en los puntos (2,16) a  (4,-5) hay un paso de una region de Y positiva a una negativa, 
por lo cual debio haber pasado por el eje X
"""

def biseccion(func, a, b, tol):
    iterac = 0
    while(b-a)/2 >= tol:
        medio = (a+b)/2
        if func(medio) == 0:
            return medio, iterac
        elif func(a) * func(medio) < 0:
            b = medio
        else:
            a = medio
        iterac += 1
    return (a + b) / 2, iterac

raiz, iterac = biseccion(Spline, 2, 4, 10e-5)

print(f"La raiz positiva es: {raiz}")

f = lambda t,y: -2*(np.sin(2*t - np.exp(-t**2 * y)))

def heun(f, t0, y0, h, n_pasos):
    """Método de Heun (Euler mejorado / RK2)"""
    t = t0
    y = y0
    t_vals = [t]
    y_vals = [y]
    for i in range(n_pasos):
        k1 = f(t, y)
        k2 = f(t + h, y + h * k1)
        y = y + (h / 2) * (k1 + k2)
        t = t + h
        t_vals.append(t)
        y_vals.append(y)
        
    return np.array(t_vals), np.array(y_vals)

inter = [0,3]

t0 = 0
y0 = 1
h = (inter[1] - inter[0])/30
n = 30

t_h,y_h = heun(f,t0,y0,h,n)

integral_Trapecios1 = np.trapezoid(y_h, t_h)
print(f"Integral con Trapecios: {integral_Trapecios1:.6f}")

def simpson(x, y):
    n = len(x) - 1 
    h = (x[len(x) - 1] - x[0]) / n
    
    if n % 2 == 0:  # n debe ser PAR para Simpson
        suma_impares = sum(y[i] for i in range(1, n, 2))
        suma_pares = sum(y[i] for i in range(2, n-1, 2))
        integral = (h/3) * (y[0] + y[n] + 4*suma_impares + 2*suma_pares)
    else:
        # Si n es impar, usar trapecio
        print(f"""n es impar ({n}), aplicando método de Trapecios para el último segmento.""")
        integral = simpson(x[:len(x) - 1], y[:len(y) - 1]) + (h/2) * (y[len(y) - 2] + y[len(y) - 1])
    
    return integral

integral_Simpson = simpson(t_h, y_h)
print(f"Integral con Simpson: {integral_Simpson:.6f}")

print("si, se puede utilizar Simpson dado que el numero de pasos es par")

h_200 = (inter[1] - inter[0])/200
n_200 = 200

t_h_200,y_h_200 = heun(f,t0,y0,h_200,n_200)

integral_Trapecios2 = np.trapezoid(y_h_200, t_h_200)
print(f"Integral con Trapecios 200 pasos: {integral_Trapecios2:.6f}")
integral_Simpson2 = simpson(t_h_200, y_h_200)
print(f"Integral con Simpson 200 pasos: {integral_Simpson2:.6f}")

plt.axhline(0, color='k', linewidth=1)
plt.plot(t_h, y_h, label='Heun 30 pasos')
plt.plot(t_h_200, y_h_200, label='Heun 200 pasos')
plt.xlabel("T")
plt.ylabel("Y")
plt.legend()
plt.grid(True)
plt.show()



