##  Bienvenido a ODEs Documentation

Esta pagina contiene el proyecto de documentación de las funciones para resolver ODE: Euler, RK2 y RK4.

##Table of Contents
Hecho por Manfred Alvarado López C10318

1. [Guía de Uso](reference.md)

## Vista del Proyecto


::: ODE


## Método de Euler

El método de Euler se basa en utilizar la siguiente expansión de Taylor:


$$
\\ x(t+h) = x(t) + h\frac{dx}{dt} + \overbrace{ \frac{h^2}{2} \frac{d^2x}{dt^2} } ^{\epsilon} + O(h^3)
$$

h siendo un cambio lo suficientemente pequeño como para asegurar que para avanzar en el intervalo de variables independientes, se puede utilizar la ecuación: 

$$
\\ x(t + h) = x(t) + hf(x,t)
$$

Cabe destacar que el valor de $$ \\epsilon $$ de la expansión de Taylor se relaciona al error. Este está relacionado de manera directa con la cantidad de subdivisiones que se hagan del intervalo de las variables independientes. El calculo de esta subdivisiones se relaciona a la expresión $$ \\ N = \frac{(b - a)}{h} $$ por lo que entre más subdivisiones, el error será más despreciable.


## Método de Runge-Kutta de Segundo Orden (RK2)


