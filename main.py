# Importando las bibliotecas necesarias
from typing import Tuple
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.lines import Line2D
from tqdm import trange

# Constantes
G = 6.673e-11                 # Constante gravitacional
AU = 1.496e11                 # Unidad astronómica en km
YEAR = 365*24*60*60.0         # Segundos en un año
MM = 6e24                     # Normalizando masa
ME = 6e24/MM                  # Masa normalizada de la Tierra
MS = 2e30/MM                  # Masa normalizada del Sol
MJ = 500*1.9e27/MM            # Masa normalizada de Júpiter
GG = (MM*G*YEAR**2)/(AU**3)   # Constante gravitacional para la simulación

def gravitational_force(m1: float, m2: float, r: np.ndarray) -> np.ndarray:
    F_mag = GG * m1 * m2 / (np.linalg.norm(r) + 1e-20)**2
    theta = np.arctan2(np.abs(r[1]), np.abs(r[0]) + 1e-20)
    F = F_mag * np.array([np.cos(theta), np.sin(theta)])
    F *= -np.sign(r)
    return F

# Resolutor RK4
def RK4Solver(t: float, r: np.ndarray, v: np.ndarray, h: float, planet: str, r_other: np.ndarray, v_other: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Resolutor de cuarto orden de Runge-Kutta para el movimiento planetario.
    """
    def dr_dt(v: np.ndarray) -> np.ndarray:
        return v

    def dv_dt(r: np.ndarray, planet: str) -> np.ndarray:
        if planet == 'earth':
            return (gravitational_force(ME, MS, r) + gravitational_force(ME, MJ, r - r_other)) / ME
        elif planet == 'jupiter':
            return (gravitational_force(MJ, MS, r) - gravitational_force(MJ, ME, r - r_other)) / MJ

    k11 = dr_dt(v)
    k21 = dv_dt(r, planet)

    k12 = dr_dt(v + 0.5 * h * k21)
    k22 = dv_dt(r + 0.5 * h * k11, planet)

    k13 = dr_dt(v + 0.5 * h * k22)
    k23 = dv_dt(r + 0.5 * h * k12, planet)

    k14 = dr_dt(v + h * k23)
    k24 = dv_dt(r + h * k13, planet)

    y0 = r + h * (k11 + 2 * k12 + 2 * k13 + k14) / 6
    y1 = v + h * (k21 + 2 * k22 + 2 * k23 + k24) / 6

    return y0, y1

# Configuración de la animación
def setup_animation() -> Tuple[py.Figure, py.Axes, Line2D, Line2D, py.Text]:
    """
    Configura el trazado de la animación.
    """
    # Creando una Figura y Ejes de Trazado
    fig, ax = py.subplots()

    # Estableciendo los Límites y las Marcas de los Ejes
    ax.axis('square')
    ax.set_xlim((-7.2, 7.2))
    ax.set_ylim((-7.2, 7.2))
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    # Trazando el Sol
    ax.plot(0, 0, 'o', markersize=9, markerfacecolor="#FDB813",
            markeredgecolor="#FD7813")

    # Inicializando Líneas para la Tierra y Júpiter
    line_earth, = ax.plot([], [], 'o-', color='#d2eeff',
                          markevery=10000, markerfacecolor='#0077BE', lw=2)
    line_jupiter, = ax.plot([], [], 'o-', color='#e3dccb', markersize=8,
                            markerfacecolor='#f66338', lw=2, markevery=10000)
    # Agregando un Objeto de Texto
    ttl = ax.text(0.24, 1.05, '', transform=ax.transAxes, va='center')

    # Devolviendo los Componentes
    return fig, ax, line_earth, line_jupiter, ttl

# Función de animación
def animate(i: int) -> Tuple[Line2D, Line2D, py.Text]:
    """
    Función de animación para el movimiento planetario.
    """
    earth_trail, jupiter_trail = 40, 200
    tm_yr = 'Tiempo transcurrido ='+ str(round(t[i], 1)) +'años'
    ttl.set_text(tm_yr)
    line_earth.set_data(r[i:max(1, i - earth_trail):-1, 0],
                        r[i:max(1, i - earth_trail):-1, 1])
    line_jupiter.set_data(r_jupiter[i:max(
        1, i - jupiter_trail):-1, 0], r_jupiter[i:max(1, i - jupiter_trail):-1, 1])
    return line_earth, line_jupiter, ttl

# Inicialización
ti, tf = 0, 120  # Tiempo inicial y final en años
N = 100 * tf     # 100 puntos por año
t = np.linspace(ti, tf, N)  # Arreglo de tiempo
h = t[1] - t[0]  # Paso de tiempo

# Inicialización de posición y velocidad
r = np.zeros([N, 2])         # Posición de la Tierra
v = np.zeros([N, 2])         # Velocidad de la Tierra
r_jupiter = np.zeros([N, 2])  # Posición de Júpiter
v_jupiter = np.zeros([N, 2])  # Velocidad de Júpiter

# Condiciones iniciales
r[0] = [1496e8 / AU, 0]
r_jupiter[0] = [5.2, 0]
v[0] = [0, np.sqrt(MS * GG / r[0, 0])]
v_jupiter[0] = [0, 13.06e3 * YEAR / AU]


# Ejecutando la simulación
for i in trange(N - 1, desc="Generando Animación"):
    r[i + 1], v[i + 1] = RK4Solver(t[i], r[i],
                                   v[i], h, 'earth', r_jupiter[i], v_jupiter[i])
    r_jupiter[i + 1], v_jupiter[i +
                                1] = RK4Solver(t[i], r_jupiter[i], v_jupiter[i], h, 'jupiter', r[i], v[i])


# Configurando la animación
fig, ax, line_earth, line_jupiter, ttl = setup_animation()
# Agregando escala y etiquetas
ax.plot([-6,-5],[6.5,6.5],'r-')
ax.text(-4.5,6.3,r'1 UA = $1.496 \times 10^8$ km')

ax.plot(-6,-6.2,'o', color = '#d2eeff', markerfacecolor = '#0077BE')
ax.text(-5.5,-6.4,'Tierra')

ax.plot(-3.3,-6.2,'o', color = '#e3dccb',markersize = 8, markerfacecolor = '#f66338')
ax.text(-2.9,-6.4,'Super Júpiter (500x masa)')

ax.plot(5,-6.2,'o', markersize = 9, markerfacecolor = "#FDB813",markeredgecolor ="#FD7813")
ax.text(5.5,-6.4,'Sol')

# Creando la animación
anim = animation.FuncAnimation(
    fig, animate, frames=4000, interval=1, blit=False)

# Mostrando la animación
plt.show()