"""
Ejercicio 1: Dada la función f (x) = ex4 4x3+2x2+8x 8 +x 3 Consideremos la ecuación x 2 f (x) = 1 
1. Haga un dibujo en el intervalo 0; 9/5 de una función convenientemente definida, de manera tal que sus raíces sean las soluciones de la ecuación dada.
¿Cuántas soluciones hay en ese intervalo? (puede responder esto en base al gráfico hecho). 
2. ¿Se puede aplicar el método de bisección?. En tal caso, calcular la solución de la ecuación con 4 dígitos decimales exactos.
¿Cuántos pasos de bisección fueron necesarios para obtener esa exactitud? 
3. ¿Se puede aplicar Newton?, ¿conviene hacerlo?. Explique cuál es la ventaja o la desventaja de aplicar Newton.
"""

#Bibliotecas
import numpy as np
import matplotlib.pyplot as plt

#Ecuac
f = lambda x: (np.exp(x**4 - 4*x**3 + 2*x**2 + 8*x - 8) + x - 3)/(x-2)
g = lambda x: f(x) - 1
gprima = lambda x: ((np.exp(-4*x**3 - 8)*((4*x**4 - 20*x**3 + 28*x**2 - 17)*np.exp(x**4 + 2*x**2 + 8*x) + np.exp(4*x**3 + 8)))/((x - 2)**2))

#Variables
inter = [0,9/5] 
tol = 1e-4
x = np.linspace(inter[0], inter[1], 1000)
y = g(x)

#Grafico
plt.axhline(0, color='k', linewidth=1)
plt.plot(x, y, label='g(x) = f(x) - 1')
plt.legend()
plt.grid(True)
plt.show()

"""
# Respuesta 1:
en el intervalo se observa una única solución, ya que la función g(x) cruza el eje x una sola vez.
"""

#Funciones
def bolzano(func, a, b):
    """Verifica el cumplimiento del teorema de Bolzano."""
    return func(a) * func(b) < 0

def biseccion(func, a, b, tol):
    iterac = 0
    while(b-a)/2 > tol:
        medio = (a+b)/2
        if func(medio) == 0:
            return medio, iterac
        elif func(a) * func(medio) < 0:
            b = medio
        else:
            a = medio
        iterac += 1
    return (a + b) / 2, iterac

#Implementacion
if bolzano(g,inter[0],inter[1]):
    resultado, iterac = biseccion(g,inter[0], inter[1], tol)
    print("La raiz por biseccion da: ",round(resultado,4),"y tuvo iteraciones", iterac)
else:
    print("El intervalo no sirve para usar biseccion")


"""
1.2. Si, se puede aplicar el metodo de biseccion. se requieren 14 iteraciones para llegar a 4 digitos de exactitud
"""

def newtonraphson(func,dfunc, a, tol):
    x0 = a
    iterac = 0
    while True:
        x1 = x0 - func(x0) / dfunc(x0)
        iterac += 1
        if abs(x1 - x0) < tol:
            break
        x0 = x1
    return x1, iterac

resultado,iterac=newtonraphson(g,gprima,1.2,tol)
print("La raiz por NR da: ",round(resultado,4),"y tuvo iteraciones", iterac)

#Fin del ejercicio 1

#Funciones
f2 = lambda x: np.sin(np.log(3*x))
f2prima = lambda x: (np.cos(np.log(3*x)))/x

#Variables
x = np.linspace(1, 10, 1000)
y = f2(x)

#Grafico
plt.axhline(0, color='k', linewidth=1)
plt.plot(x, y, label='f(x)')
plt.legend()
plt.grid(True)
plt.show()
raiz,iterac=newtonraphson(f2,f2prima,3,tol)
print("La raiz por NR da: ",round(raiz,4),"y tuvo iteraciones", iterac)

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

inter = [1, raiz]
nodos = 7
x = np.linspace(inter[0], inter[1], nodos)
y = f2(x)

integral_Trapecios = np.trapezoid(y, x)
print(f"Integral con Trapecios: {integral_Trapecios:.6f}")
integral_simpson = simpson(x, y)
print(f"Integral con Simpson: {integral_simpson:.6f}")

fprima2 = lambda x: -( (np.sin(np.log(3*x)) + np.cos(np.log(3*x)) ) / (x**2) )
fprima4 = lambda x: -((10*np.sin(np.log(3*x)))/(x**4))

# Calcular directamente

sub = nodos - 1
h = (inter[1] - inter[0]) / sub

x_eval = np.linspace(inter[0], inter[1], 1000)
M2 = np.max(np.abs(fprima2(x_eval)))
M4 = np.max(np.abs(fprima4(x_eval)))

error_trap = ((inter[1] - inter[0]) * h**2 / 12) * M2
error_simp = ((inter[1] - inter[0]) * h**4 / 180) * M4

print(f"Cota error Trapecio: {error_trap:.10f}")
print(f"Cota error Simpson: {error_simp:.10f}")

#Fin del ejercicio 2

from scipy.interpolate import CubicSpline

X = np.array([0,1,2,3,5,7])
Y = np.array([1,3,8,7,2,7])
nodos = 8

Spline = CubicSpline(X, Y)
x_fino = np.linspace(min(X), max(X), nodos)
y_spline = Spline(x_fino)

#Grafico
plt.axhline(0, color='k', linewidth=1)
plt.plot(X, Y, 'r*', markersize=8, label='Datos')
plt.plot(x_fino, y_spline, 'g-', label='Spline cúbico')
plt.legend()
plt.grid(True)
plt.show()

integral_Trapecios = np.trapezoid(y_spline, x_fino)
print(f"Integral con Trapecios: {integral_Trapecios:.6f}")
integral_simpson = simpson(x_fino, y_spline)
print(f"Integral con Simpson: {integral_simpson:.6f}")

inter = [0,7]
sub = nodos - 1

# Primera derivada
dy_dx = np.gradient(y_spline, x_fino)

# Segunda derivada (derivada de la primera derivada)
d2y_dx2 = np.gradient(dy_dx, x_fino)

x_eval = np.linspace(inter[0], inter[1], 1000)
M2 = np.max(np.abs(d2y_dx2))

h = (inter[1] - inter[0]) / sub
error_trap = ((inter[1] - inter[0]) * h**2 / 12) * M2
print(f"Cota error Trapecio: {error_trap:.10f}")

#Fin del ejercicio 3

f = lambda t,y: np.exp(np.sin(t + y)) - 1

# Verificar punto de equilibrio en (π, 0)
t_equilibrio = np.pi
y_equilibrio = 0
valor_derivada = f(t_equilibrio, y_equilibrio)
print(f"f(π, 0) = {valor_derivada:.10f}")
print(f"¿Es punto de equilibrio? {np.abs(valor_derivada) < 1e-10}")

# La solución con y(π)=0 es constante: y(t) = 0 para todo t
print("Solución con y(π)=0: y(t) = 0 para todo t")

def euler(f, t0, y0, h, n_pasos):
    """Método de Euler explícito"""
    t = t0
    y = y0
    t_vals = [t]
    y_vals = [y]
    
    for i in range(n_pasos):
        y = y + h * f(t, y)
        t = t + h
        t_vals.append(t)
        y_vals.append(y)
    
    return np.array(t_vals), np.array(y_vals)

def heun(f, t0, y0, h, n_pasos):
    """Método de Euler modificado (Heun)"""
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

# Parámetros
t0 = np.pi
y0_nuevo = 0.1
t_final = 2 * np.pi

# a) Euler con 1 paso
h_euler_1 = t_final - t0
t_euler_1, y_euler_1 = euler(f, t0, y0_nuevo, h_euler_1, 1)
print(f"Euler 1 paso: y(2π) = {y_euler_1[-1]:.6f}")

# b) Euler con 4 pasos
n_euler_4 = 4
h_euler_4 = (t_final - t0) / n_euler_4
t_euler_4, y_euler_4 = euler(f, t0, y0_nuevo, h_euler_4, n_euler_4)
print(f"Euler 4 pasos: y(2π) = {y_euler_4[-1]:.6f}")

# c) Heun con 1 paso
h_heun_1 = t_final - t0
t_heun_1, y_heun_1 = heun(f, t0, y0_nuevo, h_heun_1, 1)
print(f"Heun 1 paso: y(2π) = {y_heun_1[-1]:.6f}")

# d) Heun con 4 pasos
n_heun_4 = 4
h_heun_4 = (t_final - t0) / n_heun_4
t_heun_4, y_heun_4 = heun(f, t0, y0_nuevo, h_heun_4, n_heun_4)
print(f"Heun 4 pasos: y(2π) = {y_heun_4[-1]:.6f}")

