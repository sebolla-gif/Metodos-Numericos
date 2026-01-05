import numpy as np
import matplotlib.pyplot as plt


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

def euler_modificado(f, t0, y0, h, n_pasos):
    """Euler modificado (RK2) - método del punto medio"""
    t = t0
    y = y0
    t_vals = [t]
    y_vals = [y]

    for i in range(n_pasos):
        k1 = f(t, y)
        k2 = f(t + h/2, y + (h/2)*k1)
        y = y + h * k2
        t = t + h
        t_vals.append(t)
        y_vals.append(y)

    return np.array(t_vals), np.array(y_vals)

def rk4(f, t0, y0, h, n_pasos):
    t = t0
    y = y0
    t_vals = [t]
    y_vals = [y]

    for i in range(n_pasos):
        k1 = f(t, y)
        k2 = f(t + h/2, y + h*k1/2)
        k3 = f(t + h/2, y + h*k2/2)
        k4 = f(t + h,   y + h*k3)
        y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        t = t + h
        t_vals.append(t)
        y_vals.append(y)

    return np.array(t_vals), np.array(y_vals)

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

t0, y0 = 0, 0.5
tf = 10
n = 200   # podes cambiarlo para ver diferencias
h = (tf - t0) / n

# Ejecutar los métodos
t_euler, y_euler = euler(f, t0, y0, h, n)
t_heun, y_heun = heun(f, t0, y0, h, n)
t_emod, y_emod = euler_modificado(f, t0, y0, h, n)
t_rk4, y_rk4 = rk4(f, t0, y0, h, n)

# Gráfico comparativo
plt.figure(figsize=(10,6))
plt.plot(t_rk4, y_rk4, label="RK4 (referencia)", linewidth=3)
plt.plot(t_heun, y_heun, label="Heun (RK2)")
plt.plot(t_emod, y_emod, label="Euler Modificado (midpoint)")
plt.plot(t_euler, y_euler, label="Euler (RK1)", linestyle="--")

plt.grid(True)
plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Comparación de Euler – Euler Modificado – Heun – RK4")
plt.legend()
plt.show()
