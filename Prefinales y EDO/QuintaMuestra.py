#Bibliotecas
import numpy as np
import matplotlib.pyplot as plt

#Ecuac
f = lambda x: ((np.sqrt((x**2) + 1)) + ((1/3) * x**3) - 1) -x
fprima = lambda x: x/(np.sqrt(x**2 + 1)) + x**2 -1

# x=0 cumple con lo pedido

#Variables
inter = [-1,2] 
tol = 1e-6
x = np.linspace(inter[0], inter[1], 10000)
y = f(x)

#Grafico
plt.axhline(0, color='k', linewidth=1)
plt.plot(x, y, label='g(x) = f(x) - x')
plt.legend()
plt.grid(True)
plt.show()

#Funciones
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

resultado,iterac=newtonraphson(f,fprima,1,tol)
print("La raiz por NR da: ",resultado,"y tuvo iteraciones", iterac)

g = lambda x: np.cos(np.sqrt( 1 + x**2 ))
gprima = lambda x: -(x*np.sin(np.sqrt(x**2 + 1)))/(np.sqrt(x**2 + 1))

#Variables
inter = [0,6] 
tol = 1e-6
x = np.linspace(inter[0], inter[1], 1000)
y = g(x)

#Grafico
plt.axhline(0, color='k', linewidth=1)
plt.plot(x, y, label='g(x)')
plt.legend()
plt.grid(True)
plt.show()

raiz1,iterac=newtonraphson(g,gprima,1,tol)
print("La raiz por NR da: ",raiz1,"y tuvo iteraciones", iterac)
raiz2,iterac=newtonraphson(g,gprima,5,tol)
print("La raiz por NR da: ",raiz2,"y tuvo iteraciones", iterac)

x = np.linspace(inter[0], inter[1], 5000)

y = gprima(x)

g2prima = np.gradient(y, x)

# Valor absoluto para la cota
abs_fpp_numeric = np.abs(g2prima)
M2 = np.max(abs_fpp_numeric)

print(f"Cota numérica de |f''(x)|: M₂ = {M2:.6f}")

#Grafico
plt.axhline(0, color='k', linewidth=1)
plt.plot(x, abs_fpp_numeric, label='g(x)')
plt.legend()
plt.grid(True)
plt.show()

nodos = 2

inter = [(raiz1-0),(6-raiz2)] 

longtotal = inter[1] + inter[0]

while(True):

    sub = nodos - 1
    h = (longtotal) / sub
    error_trap = ((longtotal) * h**2 / 12) * M2
    if error_trap < 1e-4:
        print(f"Cota error Trapecio: {error_trap:.10f}")
        print(nodos)
        break
    nodos += 1

inter = [0,raiz1]

nodosinterv1 = nodos/((longtotal)/(inter[1] - inter[0]))

x = np.linspace(inter[0], inter[1], int(nodosinterv1))
y = g(x)
integral_Trapecios1 = np.trapezoid(y, x)
print(f"Integral con Trapecios: {integral_Trapecios1:.6f}")

inter = [raiz2,6]
nodosinterv2 = nodos/((longtotal)/(inter[1] - inter[0]))
x = np.linspace(inter[0], inter[1], int(nodosinterv2))
y = g(x)
integral_Trapecios2 = np.trapezoid(y, x)
print(f"Integral con Trapecios: {integral_Trapecios2:.6f}")

areatotal = integral_Trapecios1 + integral_Trapecios2
print(f"Integral con Trapecios total: {areatotal:.6f}")

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

f = lambda t,y: np.log(t + y**2) - y
t0= 0
y0 = 0.5
tf = 10
h = tf - t0

#10 n
t_euler_1, y_euler_1 = euler(f, t0, y0, h/10, 10)
print(f"Euler 10 paso: y(2π) = {y_euler_1[-1]:.6f}")

plt.axhline(0, color='k', linewidth=1)
plt.plot(t_euler_1, y_euler_1, label='Euler 10 pasos')
plt.legend()
plt.grid(True)
plt.show()

#20 n
t_euler_2, y_euler_2 = euler(f, t0, y0, h/20, 20)
print(f"Euler 20 paso: y(2π) = {y_euler_2[-1]:.6f}")

plt.axhline(0, color='k', linewidth=1)
plt.plot(t_euler_2, y_euler_2, label='Euler 20 pasos')
plt.plot(t_euler_1, y_euler_1, label='Euler 10 pasos')
plt.legend()
plt.grid(True)
plt.show()

#50 n
t_euler_5, y_euler_5 = euler(f, t0, y0, h/50, 50)
print(f"Euler 50 paso: y(2π) = {y_euler_5[-1]:.6f}")

plt.axhline(0, color='k', linewidth=1)
plt.plot(t_euler_5, y_euler_5, label='Euler 50 pasos')
plt.plot(t_euler_2, y_euler_2, label='Euler 20 pasos')
plt.plot(t_euler_1, y_euler_1, label='Euler 10 pasos')
plt.legend()
plt.grid(True)
plt.show()

#100 n
t_euler_10, y_euler_10 = euler(f, t0, y0, h/100, 100)
print(f"Euler 100 paso: y(2π) = {y_euler_10[-1]:.6f}")

plt.axhline(0, color='k', linewidth=1)
plt.plot(t_euler_10, y_euler_10, label='Euler 100 pasos')
plt.plot(t_euler_5, y_euler_5, label='Euler 50 pasos')
plt.plot(t_euler_2, y_euler_2, label='Euler 20 pasos')
plt.plot(t_euler_1, y_euler_1, label='Euler 10 pasos')
plt.legend()
plt.grid(True)
plt.show()

ft = lambda t,y: 1/(t+y**2)
fy = lambda t,y: (2*y/(y**2 + t)) - 1
#f2prima = ft + fy * f
# Usar la solución de Euler con n=100 como referencia
t_vals, y_vals = euler(f, 0, 0.5, h/100, 100)  # h=0.1, n=100

# Evaluar las funciones a lo largo de la solución
f_vals = np.abs(f(t_vals, y_vals))
ft_vals = np.abs(ft(t_vals, y_vals))  # ft = 1/(t + y²)
fy_vals = np.abs(fy(t_vals, y_vals))  # fy = (2y)/(t + y²) - 1

# Graficar
plt.axhline(0, color='k', linewidth=1)
plt.plot(t_vals, f_vals, label='|f(t,y)|')
plt.plot(t_vals, ft_vals, label='|f_t(t,y)|') 
plt.plot(t_vals, fy_vals, label='|f_y(t,y)|')
plt.legend()
plt.grid(True)
plt.show()

# Para acotar y''(t) necesitamos: y'' = f_t + f_y * f
def y_segunda(t, y):
    return ft(t, y) + fy(t, y) * f(t, y)

n=100
t_euler_n, y_euler_n = euler(f, t0, y0, h/n, n)
# Calculamos y'' en estos puntos
y_segunda_vals = y_segunda(t_euler_n, y_euler_n)
M = np.max(np.abs(y_segunda_vals))
tol = 1e-4
n_global = int( np.ceil((M*10) / (2*tol)) )
print("n necesario según la cota global =", n_global)


"""
Ese es el error local, pero NO sirve para decidir n, porque lo que querés controlar es el error global, que para Euler es:

Eglobal ≤ Mh/2 donde M = max |y''(t)| = max |ft + fy*f|
Entonces la cota correcta para el error global es:

Mh/2 ≤ 10e-4 ⇒ h ≤ 2e-4/M

y como:

h=10/n

entonces:

n≥10M/2⋅10e-4

code: 

def euler(f, t0, y0, h, n_pasos):
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


# Definiciones
f = lambda t,y: np.log(t + y**2) - y
ft = lambda t,y: 1/(t+y**2)
fy = lambda t,y: (2*y)/(t+y**2) - 1
y2 = lambda t,y: ft(t,y) + fy(t,y)*f(t,y)

t0=0
y0=0.5
tf=10
h_total = tf - t0


# ---- Usar la solución con 100 pasos para estimar M ----
t_vals, y_vals = euler(f, t0, y0, h_total/100, 100)

M = np.max(np.abs(y2(t_vals, y_vals)))  # Máximo de |y''|
print("M estimado =", M)

# ---- Cálculo correcto del n usando la cota del error global de Euler ----
tol = 1e-4

# Cota global:  (M*h)/2 <= tol
# h = 10/n
# (M*(10/n))/2 <= tol
# n >= (M*10) / (2*tol)

n_global = int(np.ceil( (M*10) / (2*tol) ))
print("n necesario según la cota global =", n_global)

# ---- Gráfico de f, ft, fy ----
plt.axhline(0, color='k', linewidth=1)
plt.plot(t_vals, np.abs(f(t_vals, y_vals)), label='|f|')
plt.plot(t_vals, np.abs(ft(t_vals, y_vals)), label='|f_t|')
plt.plot(t_vals, np.abs(fy(t_vals, y_vals)), label='|f_y|')
plt.legend()
plt.grid()
plt.show()

"""

