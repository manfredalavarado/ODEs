##  Bienvenido a ODEs Documentation

Esta pagina contiene el proyecto de documentación de las funciones para resolver ODE: Euler, RK2 y RK4.

##Tabla de Contenidos

1. [Guía de Uso](reference.md)

## Vista del Proyecto


::: ODE


## Método de Euler

El método de Euler se basa en utilizar la siguiente expansión de Taylor:


$$
\\ x(t+h) = x(t) + h\frac{dx}{dt} + \overbrace{ \frac{h^2}{2} \frac{d^2x}{dt^2} } ^{\epsilon} + O(h^3)
$$


Donde h es un número de cambio lo suficientemente pequeño como para estimar la solución de la EDO de la siguiente forma:

$$
\\ x(t + h) = x(t) + h\frac{dx}{dt} =  x(t) + hf(x,t)
$$

Donde el error global escala linealmente con h:

$$
\\ \sum\epsilon = \frac{h}{2}\left[f_b - f_a\right]
$$


## Método de Runge-Kutta de Segundo Orden (RK2)

En este método, se llega a una aproximación a la solución de la forma:

$$
\\ x(t + h) = x(t) + hf[x(t + h/2), t + h/2] + O(h^3)
$$

El procedimiento usual para resolver una EDO por este método consiste en calcular los incrementos:

$$
\\ k_1 = hf(x,t),  \quad \quad \quad k_2 = hf(x+\frac{k_1}{2}, t + \frac{h}{2})
$$

Para calcular la solución evaluada en un tiempo t + h con la siguiente ecuación:

$$
\\ x(t+h) = x(t) + k_2
$$

Donde se obtiene un error global que escala con h al cuadrado.


## Método de Runge-Kutta de Cuarto Orden (RK4)

Si se busca una mayor precisión se puede utilizar el método RK4, el cual suele ser método predilecto para estos problemas.

El procedimiento consiste, al igual que con RK2, en encontrar los incrementos, pero esta vez son 4:

$$
\\ k_1 = hf(x, t), \quad \quad \quad k_2 = hf\left(x + \frac{k_1}{2}, t+\frac{h}2\right),
$$

$$
k_3 = hf\left(x + \frac{k_2}{2}, t+\frac{h}2\right), \quad \quad \quad k_4 = hf\left(x + k_3, t + h \right)
$$


Para finalmente calcular la solución evaluada en un tiempo t + h con la siguiente ecuación:

$$
\\ x(t+h) = x(t) + \frac{1}{6}(k_1 + 2 k_2 + 2k_3 + k_4).
$$
