---
title: "SIMP método para optimización topológica"
date: 2023-03-20T11:48:52-05:00
draft: false
math: true
---

Este método ha demostrado ser utilizado en una gran cantidad de problemas, ya que es un método eficiente computacionalmente y simple de implementar. A diferencia de BESO [ESO/BESO](https://kssgarcia.github.io/es/posts/eso_beso/) este método utiliza un criterio *soft-kill*, es decir que no se eliminan completamente los elementos, solo se reduce el modulo elástico o una de las características dimensionales. Este criterio se implementa ya que la eliminación completa de elementos puede resultar en dificultades teóricas.

En el problema que se presenta se busca minimizar el cumplimiento mínimo del modelo. Este esta dado por

  
$$
\begin{equation}
\text{Minimizar } C = f^T u
\tag{1a}
\end{equation}
$$

$$
\begin{equation}
\text{Sujeto a: } V^* - \sum _{i=1}^N v_i x_i = 0
\tag{1b}
\end{equation}
$$

$$
\begin{equation}
x_i = x_{min} \text{ o } 1
\tag{1c}
\end{equation}
$$


$v_i$ es el volumen del elemento, $V^*$ es el volumen elegido para la estructura final y $x_i$ es la variable de diseño que denota la densidad relativa del elemento $i$; en este método se utiliza el doble de la energía de deformación que se utiliza en el método [ESO/BESO](https://kssgarcia.github.io/es/posts/eso_beso/). Al contrario que el método BESO, la variable de diseño nunca es cero, por lo tanto nunca se eliminan los elementos completamente, evitando tener elementos vacíos; para los elementos con $x_i=x_{min}$ se les llama elementos suaves o *soft-elements* en ingles, el valor de la variable $x_{min}$ puede ser por ejemplo $0.001$ pero nunca cero.

En este método se utiliza un esquema de interpolación de material con penalización para dirigir la solución a diseños casi solidos. Para conseguir este tipo de diseños el modulo de elasticidad es interpolado como función de la densidad del elemento o la variable de diseño:


$$
\begin{equation}
E(x_i) = E_1x_i^p
\tag{2}
\end{equation}
$$

  
donde $E_1$ es el modulo de young del elemento del material solido y $p$ es el exponente de penalización. Se asume que la razón de Poisson es independiente de la variable de diseño y la matriz de rigidez global es expresada por la matriz del elemento y la variable de diseño $x_i$ como
  

$$
\begin{equation}
K = \sum_i x_i^p K_i^0
\tag{3}
\end{equation}
$$

  
donde $K^0_i$ es la matriz de rigidez del elemento solido.

### Análisis de sensibilidad 

Para obtener el numero de sensibilidad primero se considera la ecuación de equilibrio que esta expresada por


$$
\begin{equation}
Ku = f
\tag{4}
\end{equation}
$$


si se considera el cumplimiento mínimo (1a) y se asume que la variable de diseño $x_i$ varia continuamente entre $x_{min}$ y $1$, la sensibilidad de la función objetivo con respecto a la variable de diseño es
  

$$
\begin{equation}
\frac{dC}{dx_i} = \frac{df^T}{dx_i}u + f^T\frac{du}{dx_i}
\tag{5}
\end{equation}
$$
  

El método adjunto es utilizado para calcular la sensibilidad del vector de desplazamientos. Introduciendo un vector de multiplicadores de Lagrange $\lambda$, un termino extra $\lambda^T(f-Ku)$ puede ser añadido a la función objetivo sin cambiar nada debido a la ecuación (4)
  
$$
\begin{equation}
C = f^Tu + \lambda^T(f-Ku)
\tag{6}
\end{equation}
$$

teniendo en cuenta la variación de la variable de diseño
  

$$
\begin{equation}
\frac{dC}{dx_i} = \frac{df^T}{dx_i}u + f^T\frac{du}{dx_i} + \frac{\lambda^T}{dx_i}(f-Ku) + \lambda^T\left(\frac{f^T}{dx_i}-\frac{K}{dx_i}u - K\frac{u}{dx_i}\right)
\tag{7}
\end{equation}
$$


en la ecuación mostrada anteriormente se puede eliminar los términos numero 1 y 3 debido a que se asume que la variación de $x_i$ no tiene efecto en el vector de fuerzas, por lo tanto $\frac{df}{dx_i}$ y por que $f-Ku = 0$ debido a (4)
  
$$
\begin{equation}
\frac{dC}{dx_i} = (f^T - \lambda^T K) \frac{du}{dx_i} - \lambda^T \frac{dK}{dx_i} u
\tag{8}
\end{equation}
$$


de la ecuación (6) como $(f-Ku)$ es igual a cero, el vector $\lambda$ puede ser elegido como sea. Para eliminar el termino $\frac{du}{dx_i}$ de la ecuación anterior, $\lambda$ se define teniendo en cuenta


$$
\begin{equation}
f^T - \lambda^T K = 0
\tag{9}
\end{equation}
$$


con la ecuación anterior y la ecuación de equilibrio (4) se obtiene la siguiente solución para $\lambda$


$$
\begin{equation}
\lambda = u
\tag{10}
\end{equation}
$$


sustituyendo $\lambda$ en la ecuación (8),


$$
\begin{equation}
\frac{dC}{dx_i} = - u^T \frac{dK}{dx_i} u
\tag{11}
\end{equation}
$$


sustituyendo la ecuación (3) en la ecuación anterior se obtiene

  
$$
\begin{equation}
\frac{dC}{dx_i} = - p x_i^{p-1} u^T_i K_i^0 u_i
\tag{12}
\end{equation}
$$

  
El problema de optimización anterior se puede resolver por diferentes acercamiento como el del criterio de optimalidad. Un esquema para la actualización de la variable de diseño puede ser

$$
\begin{equation}
x_i^{K+1}= \begin{cases}\max \left(x_{\min }, x_i^K-m\right) & \text { si } x_i^K B_i^\eta \leq \max \left(x_{\min }, x_i^K-m\right) \\\
\min \left(1, x_i^K+m\right) & \text { si } \min \left(1, x_i^K+m\right) \leq x_i^K B_i^\eta \\\
x_i^K B_i^\eta & \text { de otro modo }\end{cases}
\tag{13}
\end{equation}
$$

donde $x_i^K$ denota el valor de diseño en la iteración $K$, $m$ es el movimiento limite positivo, $\eta$ es el coeficiente de amortiguamiento (normalmente es $0.5$) y $B_i$ se obtiene de la condición de optimalidad


$$
\begin{equation}
B_i = \lambda^{-1} p x_i^{p-1} u_i^T K_i^0 u_i
\tag{14}
\end{equation}
$$


donde $\lambda$ es el multiplicador de Lagrange que puede ser determinado usando el método de la bisección o Newton. Para asegurar que el diseño óptimo es independiente del mallado y libre del patrón de ajedrez, se introduce el siguiente esquema de filtro que modifica la sensibilidad de los elementos


$$
\begin{equation}
\frac{\partial c}{\partial x_i}=\frac{1}{x_i \sum_{j=1}^N H_{i j}} \sum_{j=1}^N H_{i j} x_j \frac{\partial c}{\partial x_j}
\tag{15}
\end{equation}
$$

  
donde $N$ es el numero total de elementos en la malla y $H_{ij}$ es el factor de peso dado por
  

$$
\begin{equation}
H_{i j}=r_{\min }-r_{i j}, \quad \lbrace i \in N | r_{ij} \leq r_{min} \rbrace
\tag{16}
\end{equation}
$$

donde $r_{ij}$ es la distancia entre los centros de los elementos $i$ y $j$. El peso $H_{ij}$ es cero fuera del circulo de radio $r_{min}$.


1. Discretizar el dominio e inicializar las variables de diseño.
2. Realice el análisis de elementos finitos.
3. Calcular la sensibilidad de la función objetivo con respecto a la densidad de cada elemento.
4. Realizar los criterios de optimalidad.
5. Calcular el error de la densidad.
6. Repita los pasos 2 a 5 hasta que se logre la masa $M^*$ o la convergencia del error.

## Referencias

- Yago, D., Cante, J., Lloberas-Valls, O., & Oliver, J. (2022). Métodos de optimización de topología para problemas estructurales 3D: un estudio comparativo. Archives of Computational Methods in Engineering, 29(3), 1525–1567. doi:10.1007/s11831-021-09626-2

- Huang, X. y Xie, M. (2010). Optimización de topología evolutiva de estructuras continuas: métodos y aplicaciones. John Wiley & Sons.