# Runge-Kutta

Simulación del movimiento orbital de la Tierra y un **Super Júpiter** (500× la masa de Júpiter) alrededor del Sol, integrada numéricamente con el método de **Runge–Kutta de cuarto orden (RK4)**.

## Descripción

El script resuelve las ecuaciones del movimiento en un sistema de dos cuerpos acoplados (Tierra–Júpiter) bajo gravitación newtoniana, con unidades normalizadas (UA, años). Al finalizar la integración muestra una animación con las trayectorias y el tiempo transcurrido en años.

- **Horizonte temporal:** 0–120 años  
- **Resolución:** 100 pasos por año  
- **Salida:** ventana interactiva de Matplotlib con animación

## Requisitos

- Python 3.10+ (recomendado)
- Dependencias en `requirements.txt`: `numpy`, `matplotlib`, `tqdm`

## Instalación

```bash
git clone https://github.com/mendozalz/Runge-Kutta.git
cd Runge-Kutta
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Uso

```bash
python main.py
```

La simulación puede tardar un poco: se integran miles de pasos (barra de progreso con `tqdm`) antes de abrir la animación.

## Estructura

| Archivo            | Descripción                                      |
|--------------------|--------------------------------------------------|
| `main.py`          | Fuerzas gravitacionales, RK4, animación          |
| `requirements.txt` | Dependencias de Python                         |
| `.gitignore`       | Excluye `venv/` y artefactos locales             |

## Licencia

Uso libre para aprendizaje y experimentación. Si reutilizas el código, un enlace al repositorio es bienvenido.
