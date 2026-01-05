import numpy as np

# ==============================
# METODOS (SIN CAMBIAR)
# ==============================

def heun(f, t0, y0, h, n_pasos):
    t = t0
    y = y0
    t_vals = [t]
    y_vals = [y]

    for _ in range(n_pasos):
        k1 = f(t, y)
        k2 = f(t + h, y + h*k1)
        y = y + (h/2)*(k1 + k2)
        t = t + h
        t_vals.append(t)
        y_vals.append(y)

    return np.array(t_vals), np.array(y_vals)

def euler_modificado(f, t0, y0, h, n_pasos):
    t = t0
    y = y0
    t_vals = [t]
    y_vals = [y]

    for _ in range(n_pasos):
        k1 = f(t, y)
        k2 = f(t + h/2, y + (h/2)*k1)
        y = y + h*k2
        t = t + h
        t_vals.append(t)
        y_vals.append(y)

    return np.array(t_vals), np.array(y_vals)

def euler(f, t0, y0, h, n_pasos):
    t = t0
    y = y0
    t_vals = [t]
    y_vals = [y]

    for _ in range(n_pasos):
        y = y + h*f(t, y)
        t = t + h
        t_vals.append(t)
        y_vals.append(y)

    return np.array(t_vals), np.array(y_vals)


# ==============================
# EDO DEL PARCIAL
# ==============================

def f(t, y):
    return np.log(t + y*y) - y

def ft(t, y):
    return 1/(t + y*y)

def fy(t, y):
    return (2*y/(t + y*y)) - 1

def y_segunda(t, y):
    return ft(t, y) + fy(t, y)*f(t, y)   # visto en la teoría


# ==============================
# PARÁMETROS
# ==============================

t0 = 0
y0 = 0.5
T = 10
N = 100
h = (T - t0)/N


# ==============================
# CÁLCULOS
# ==============================

t_eu, y_eu = euler(f, t0, y0, h, N)
t_em, y_em = euler_modificado(f, t0, y0, h, N)
t_he, y_he = heun(f, t0, y0, h, N)

# M = max |y''|
M_eu = np.max(np.abs(y_segunda(t_eu, y_eu)))
M_em = np.max(np.abs(y_segunda(t_em, y_em)))
M_he = np.max(np.abs(y_segunda(t_he, y_he)))

# ===========================
# COTAS DEL ERROR LOCAL
# ===========================

ELT_eu = (h**2)/2 * M_eu        # Euler: orden 1
ELT_em = (h**3)/6 * M_em        # EM: orden 2
ELT_he = (h**3)/6 * M_he        # Heun: orden 2

# ===========================
# COTAS DEL ERROR GLOBAL
# (SEGÚN LA TEORÍA DEL APUNTE)
# ===========================

EG_eu = M_eu * h                # O(h)
EG_em = M_em * h*h             # O(h^2)
EG_he = M_he * h*h             # O(h^2)

# ===========================
# MOSTRAR
# ===========================

print("----- COTAS DEL ERROR LOCAL -----")
print("Euler local          =", ELT_eu)
print("Euler modificado     =", ELT_em)
print("Heun                 =", ELT_he)

print("\n----- COTAS DEL ERROR GLOBAL (sin Lipschitz) -----")
print("Euler global         =", EG_eu)
print("Euler modificado     =", EG_em)
print("Heun                 =", EG_he)

print("\nNota:")
print("Estas cotas siguen EXACTAMENTE las fórmulas del tema 6 (EDOs) del apunte.")
print("La constante de Lipschitz NO se usa porque NO aparece en la teoría provista.")
