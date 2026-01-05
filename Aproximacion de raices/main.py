"""
Sebastian Vaccaro - 2do cuatri 2025
Trabajo Práctico 1: Ecuaciones No Lineales - Métodos Numéricos - UNSAM
"""

import numpy as np
import matplotlib.pyplot as plt
from metodos import resolver, bolzano

# --- La funcion y su derivada ---
f = lambda x: (3/5) * np.exp((15/7) * np.cos(3*x)) - 2*x
fprima = lambda x: (3/5) * np.exp((15/7) * np.cos(3*x)) * ((15/7) * (-3) * np.sin(3*x)) - 2

# --- Parametros ---
N = 5          # Iteraciones 
tol = 5e-6     # Tolerancia 

# --- Primer intervalo ---
intervalo1 = np.array([-9, 1])

if bolzano(f, *intervalo1):
    intervalo_final, iters, raiz = resolver(f, fprima, N, intervalo1, tol)
    print("Intervalo inicial:", intervalo1)
    print("Intervalo final (bisección):", intervalo_final)
    print("Iteraciones N-R:", iters)
    print("Raíz encontrada:", round(raiz, 5))
else:
    print("El teorema de Bolzano no se cumple en el intervalo", intervalo1)

# --- Grafico ---
x = np.linspace(-5, 5, 400)
plt.plot(x, f(x), label="f(x)")
plt.axhline(0, color="black", linewidth=0.8)
plt.title("f(x) en [-5, 5]")
plt.grid()
plt.legend()
plt.show()

# --- Segundo intervalo ---
intervalo2 = np.array([1, 2])

if bolzano(f, *intervalo2):
    intervalo_final, iters, raiz = resolver(f, fprima, N, intervalo2, tol)
    print("\nIntervalo inicial:", intervalo2)
    print("Intervalo final (bisección):", intervalo_final)
    print("Iteraciones N-R:", iters)
    print("Raíz encontrada:", round(raiz, 5))
else:
    print("El teorema de Bolzano no se cumple en el intervalo", intervalo2)
