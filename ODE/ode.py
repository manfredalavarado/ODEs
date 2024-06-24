# ODE/ode.py

"""
Contiene multiples formas de resolver ecuaciones diferenciales oridnarias.

Este modulo permite al usuario resolver sus ODEs mediante los metodos de Euler, RK2 y RK4.

"""


def Euler(f,x0,t0,tf,N):
    """
    Metodo de Euler para resolver ODE.

    Examples:
        >>> import numpy as np
        >>> def f(x,t):
        >>>     return (-1*x)**3 + np.sin(t)
        >>> x,t = Euler(f,0,0,10,20)
        >>> print(x)
        [ 0.          0.          0.26439534  0.71189381  1.04830651  0.89488942
  0.77464602  0.52141013  0.17502349 -0.28921312 -0.80263945 -0.97897518
 -0.73458315 -0.50879969 -0.16038526  0.30726691  0.81787725  0.97386719
  0.72957635  0.49945695]
        >>> print(t)
        [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
  3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
  6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
  9.47368421 10.        ]

    Args:
        f (function): Ecuación a resolver.
        x0 (int): Valor inicial para la variable dependiente.
        t0 (int): Valor inicial para la variable independiente.
        tf (int): Valor final para la variable independiente.
        N (int): Número de particiones del rango de la variable independiente.

     Return:
        Los vectores 'x' y 't', para las variables dependientes y para las variables independientes.
            

    """
    t = np.linspace(t0,tf,N)
    x = np.zeros(t.size)
    x[0] = x0
    h = t[1] - t[0]
    
    for i in range(t.size - 1):
        x[i+1] = x[i] + h*f(x[i],t[i])
    
    return x,t


def RK2(f,x0,t0,tf,N):

    """
    Método de rakuta de dos índices para resolver ODE.

    Examples:
        >>> import numpy as np
        >>> def f(x,t):
        >>>     return -x**3 + np.sin(t)
        >>> x,t = RK2(f,0,0,10,20)
        >>> print(x)
        [ 0.          0.13691106  0.50040599  0.83221858  0.89696773  0.83638552
  0.684369    0.42791827  0.0377284  -0.46988042 -0.7896376  -0.78706438
 -0.65422559 -0.40234264 -0.0089812   0.49847514  0.79691697  0.78634182
  0.64914835  0.39298319]
        >>> print(t)
        [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
  3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
  6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
  9.47368421 10.        ]

     Args:
        f (function): Ecuación a resolver.
        x0 (int): Valor inicial de la variable dependiente.
        t0 (int): Valor inicial de la variable independiente.
        tf (int): valor final de la variable independiente.
        N (int): Número de particiones del rango de la variable independiente.

    return:
       Los vectores 'x' y 't', para las variables independientes y para las variables dependientes.
    """
    
    t = np.linspace(t0,tf,N)
    x = np.zeros(t.size)
    x[0] = x0
    h = t[1] - t[0]
    
    for i in range(t.size-1):
        k1 = h*f(x[i],t[i])
        k2 = h*f(x[i] + (k1/2), t[i] + (h/2))
        n = x[i] + (h/2)*f(x[i],t[i])
        x[i+1] = x[i] + h*f(n, t[i] + h/2)
    
    return x,t



def RK4(f,x0,t0,tf,N):
    """
    Método de Rakuta con cuatro índices para resolver ODE.

    Examples:
        >>> import numpy as np
        >>> def f(x,t):
        >>>     return -x**3 + np.sin(t)
        >>> x,t = RK4(f,0,0,10,20)
        >>> print(x)
         [ 0.          0.13505937  0.48487586  0.81979626  0.93102462  0.88144889
  0.72370945  0.46026479  0.06848691 -0.42773679 -0.78084333 -0.83096351
 -0.6993985  -0.43991777 -0.04506765  0.45089026  0.79057875  0.83047637
  0.69376603  0.43014634]
        >>> print(t)
         [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
  3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
  6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
  9.47368421 10.        ]
                

    Args:
        f (func): Ecuación a resolver.
        x0 (int): Valor inicial para la variable dependiente.
        t0 (int): Valor inicial para la variable independiente.
        tf (int): Valor final para la variable independiente.
        N (int): Número de particiones del rango de la variable independiente.

    return:
        Los vectores 'x' y 'n', para las variables dependientes e independientes.
     """
    t = np.linspace(t0,tf,N)
    x = np.zeros(t.size)
    x[0] = x0
    h = t[1] - t[0]
    
    for i in range(t.size - 1):
        k1 = h*f(x[i],t[i])
        k2 = h*f(x[i] + (k1/2), t[i] + (h/2))
        k3 = h*f(x[i] + (k2/2), t[i] + (h/2))
        k4 = h*f(x[i] + k3, t[i] + h)
        x[i+1] = x[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    return x,t
