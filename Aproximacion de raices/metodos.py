"""
Sebastian Vaccaro - 2do cuatri 2025
Trabajo Práctico 1: Ecuaciones No Lineales - Métodos Numéricos - UNSAM
"""

import numpy as np

def resolver(func, func_deriv, N, intervalo, tol):

    a, b = intervalo

    # --- Metodo de Biseccion (solo N iteraciones) ---
    for _ in range(N):
        m = (a + b) / 2
        if func(a) * func(m) < 0:
            b = m
        else:
            a = m
    intervalo_final = (a, b)

    # --- Metodo de Newton-Raphson ---
    x0 = (a + b) / 2
    iteraciones_NR = 0
    while True:
        x1 = x0 - func(x0) / func_deriv(x0)
        iteraciones_NR += 1
        if abs(x1 - x0) < tol:
            break
        x0 = x1

    raiz = x1
    return intervalo_final, iteraciones_NR, raiz


def bolzano(func, a, b):
    """Verifica el cumplimiento del teorema de Bolzano."""
    return func(a) * func(b) < 0
