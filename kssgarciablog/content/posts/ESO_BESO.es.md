---
title: "ESO/BESO"
date: 2023-03-15T11:48:52-05:00
draft: false
math: true
---

El método ESO (evolutionary structural optimization) y BESO (bi-directional evolutionary structural optimization) son metodos con criterio de optimización hard-kill, es decir que puede eliminar material y esto hace que el dominio se vea modificado.

El método ESO consiste en el simple concepto de ir removiendo material que es ineficiente. Esto se realiza calculando los esfuerzos en todo el modelo y eliminando los elementos que presentan valores bajos. 

## Optimizacion por nivel de esfuerzo 

Para esto se calcula el nivel de esfuerzo de cada elemento y se compara con criterio definido como criterio de rechazo

  

$$
\begin{equation}
	\frac{\sigma_e}{\sigma_m} < RR_i 
\end{equation}
$$

  
$\sigma_e$ es el esfuerzo de Von Mises del elemento, $\sigma_m$ es el esfuerzo de Von Mises del modelo y $RR_i$ es el criterio de rechazo. Esta desigualdad define que elementos son removidos después de cada análisis de finitos. Este ciclo se repite hasta que se llegue a un estado de equilibrio, es decir que no hay mas elementos por ser removidos con el criterio de rechazo utilizado. Después de esto se puede aumentar el valor de $RR_i$ para capturar elementos con mayor nivel de esfuerzo y así lograr remover mas material

$$
\begin{equation}
	RR_{i+1} = RR_i + ER
\end{equation}
$$


$ER$ es una tasa evolutiva que se le suma al criterio de rechazo. El método ESO puede ser resumido como:

1.  Discretización del modelo en elementos finitos.
2.  Realizar un análisis de elementos finitos de la estructura.
3.  Remover elementos que satisfacen (1).
4.  Aumentar el criterio de rechazo utilizando (2).
5.  Repetir los pasos 3 y 4 hasta llegar al punto deseado.


## Optimización por energía de deformación. 

Este criterio consiste en tomar la rigidez como factor para llevar a cabo la optimización. Comúnmente llamado cumplimiento medio $C$, es la medida inversa de la rigidez de toda la estructura. El cumplimiento mínimo puede definirse como la energía de deformación total

$$
\begin{equation}
C = \frac{1}{2} f^T u
\end{equation}
$$

donde $f$ y $u$ es el vector de fuerza y desplazamientos respectivamente. El equilibrio de la estructura es definido por:

$$
\begin{equation}
Ku = f \end{equation} 
$$
donde $K$ es la matriz de rigidez global Cuando un elemento es eliminado la matriz de rigidez es reescrita de la siguiente forma:

$$K^* = K - K_i$$ $$\Delta K = K^* - K = -K_i$$


$K^*$ es la matriz modificada después de remover un elemento, $K_i$ es la matriz del elemento. En la ecuación después de eliminar un elemento solo varia el vector de desplazamiento, el cual se describe

$$
\begin{equation}
\Delta u = -K^{-1} \Delta Ku
\end{equation}
$$

de las ecuaciones (3) y (5)

  

$$
\begin{equation}
\Delta C = \frac{1}{2} f^T \Delta u = - \frac{1}{2} f^T k^{-1} \Delta K u = \frac{1}{2} u^T_i K_iu_i
\end{equation}
$$

  

donde $u_i$ es el vector de desplazamiento del elemento. $\Delta C$ se entiendo como el numero de sensibilidad del cumplimiento mínimo y es denotado por $\alpha^e_i$.

  
$$
\begin{equation}
\alpha^e_i = \frac{1}{2} u^T_i K_iu_i
\end{equation}
$$


Esta ecuación expresa que el incremento en la energía de deformación o cumplimiento mínimo es igual a la energía de deformación del elemento. Para minimizar el cumplimiento mínimo, lo cual es equivalente a maximizar la matriz de rigidez, se eliminan los elementos con menores valores de $\alpha^e_i$.

ESO con optimización de rigidez se puede resumir en los siguientes pasos:

1.  Discretización del modelo en elementos finitos.
2.  Realizar un análisis de elementos finitos de la estructura.
3.  Calcular el valor de sensibilidad de cada elemento con (7).
4.  Remover los elementos con menor valor de sensibilidad de acuerdo con el ratio de rechazo.
5.  Repetir los pasos 2 al 4 hasta que el cumplimiento mínimo o máxima deformación se consigue.

El método BESO es una versión mejorada del método ESO, ya que este además de eliminar material también añade. Los elementos solidos con un valor de sensibilidad bajo son eliminados y los elementos vacíos con valor de sensibilidad altos se añaden y se convierten en elementos solidos. El valor de sensibilidad de los elementos vacíos es calculado por medio de extrapolación. Para este método el ratio de rechazo $RR$ define el numero de elementos que son eliminados en cada iteración y el ratio de inclusión $IR$ es el que define el aumento del ratio de rechazo.

Un problema que se presenta en el método ESO es la alta sensibilidad que tiene al aplicar las restricciones, siendo posible conseguir un peor diseño que el inicial. Por lo tanto en este nuevo método el problema de optimización esta mejor planteado para que de esta forma se asegure mejor convergencia

El planteamiento del problema al que apunta este método esta definido por:

$$\begin{equation}
    Minimizar\ C = \frac{1}{2} f^T u
\end{equation}$$
$$\begin{equation}
    V^* - \sum_{i=1}^N V_i x_i = 0
\end{equation}$$

  
$V^*$ es el volumen de toda la estructura y $v_i$ es el volumen del elemento individual. $x_i$ es una variable que define cuando un elemento esta presente y cuando no lo esta, 1 o 0 respectivamente.

Para este método el numero de sensibilidad es el descrito por la expresión (7). Cuando el mallado es no uniforme, el numero de sensibilidad debe tener en cuenta el volumen del elemento. De esta forma, el numero de sensibilidad puede ser reemplazado por la densidad de energía de deformación del elemento.


$$
\begin{equation}
\alpha^e_i = \left(\frac{1}{2} u^T_i K_iu_i\right)/V_i
\end{equation}
$$

  
Cuando una estructura es discretizada utilizando pocos elementos puede generarse problemas de convergencia debido a que que el numero de sensibilidad se vuelve discontinuo. Para evitar este problema se ideo un filtro que suaviza la sensibilidad por medio de nodos y elementos para evitar la dependencia de mallado.

Para llevar el filtro de suavizado primero hay que suavizar los nodos, para ello se utiliza la siguiente ecuación
  

$$
\begin{equation}
\alpha^n_j = \sum_{i=1}^M w_i\alpha^e_i
\end{equation}
$$
  
donde $M$ es el numero de elementos conectados al nodo $j$. $w_i$ es el peso del elemento $i$ y $\sum_{i=1}^M wi = 1$.
  

$$
\begin{equation}
w_i = \frac{1}{M-1} \left( 1 - \frac{r_{ij}}{\sum_{i=1}^M r_{ij}} \right)
\end{equation}
$$

  

donde $r_{ij}$ es la distancia entre el centro del elemento $i$ y el nodo $j$, de esta forma el peso es mayor cuando el la distancia entre el elemento y el nodo es menor. Después de realizar el suavizado de los nodos se realiza el suavizad de los elementos. Para realizar este esquema se considera la distancia $r_{min}$ que no cambia con el refinamiento de la malla. Este distancia define el radio del dominio $\omega_i$. El propósito de este dominio es identificar los nodos que van a influenciar la sensibilidad del elemento $i$.
  

$$
\begin{equation}
\alpha_i = \frac{\sum_{i=1}^K w(r_{ij}) \alpha_j^n}{\sum_{i=1}^K w(r_{ij})}
\end{equation}
$$

donde K es el numero total de nodos en el dominio $\omega_i$, $w(r_{ij})$ es el factor de peso definido

  

$$
\begin{equation}
w(r_{ij}) = r_{min} - r_{ij}\ (j=1,2,...,K)
\end{equation}
$$

  

de esta forma el numero de sensibilidad es suavizado y la sensibilidad de los elementos vació es obtenida. Con este método se consigue resolver el problema generado por la dependencia de la malla, sin embargo esto genera grandes oscilaciones en la función objetivo produciendo problemas de convergencia. Para esto se encontró que al promediar la sensibilidad actual con la anterior resuelve este problema, este promedio se describe por medio de
  
$$
\begin{equation}
\alpha_i = \frac{\alpha_i^k + \alpha_i^{k-1}}{2}
\end{equation}
$$

  

donde $k$ es la iteración actual.
Para remover y añadir elementos se tiene que calcular previamente el volumen objetivo de la siguiente iteración $V_{k+1}$. Debido a que la restricción de volumen $v^*$ puede ser mayor o menor que el volumen inicial, el volumen objetivo puede aumentar o disminuir. El volumen evoluciona de la siguiente forma
  

$$
\begin{equation}
V_{k+1} = V_k(1 \pm ER)\ (k=1,2,3,...)
\end{equation}
$$

  

donde $ER$ es el ratio de evolución del volumen. Una vez la restricción de volumen es cumplida, el volumen de la siguiente iteración seguir constante.

Al calcular la sensibilidad de los elementos estos son ordenados de mayor a menor. Para los elementos solidos que serán eliminados

  

$$\alpha_i \leq \alpha_{del}^{th}$$

  

para los elementos vacíos que serán añadidos
  

$$\alpha_i > \alpha_{add}^{th}$$

  

Para calcular los valores de $\alpha_{add}^{th}$ y $\alpha_{del}^{th}$ hay que seguir los siguientes pasos:

  
1.  dejar $\alpha_{add}^{th} = \alpha_{del}^{th} = \alpha_{th}$. Si se tienen 1000 elementos y con las sensibilidades ordenadas $\alpha_1 > \alpha_2 > ... > \alpha_{1000}$ y $V_{k+1}$ corresponde a un diseño con 600 elementos entonces $\alpha_{th} = \alpha_{600}$
  

2.  Calcular la razón de adición de volumen $AR$, la cual es calcula da dividiendo el numero de elementos añadidos por el numero total de elementos. Si $AR \leq AR_{max}$ entonces no realizar el paso 3. De otra forma hay que recalcular $\alpha_{add}^{th}$ y $\alpha_{del}^{th}$.

  

3.  Para calcular $\alpha_{add}^{th}$ ordenar la sensibilidad de los elementos vacíos. El numero de elementos que van a cambiar de 0 a 1 sera igual a $AR_{max}$ multiplicado el numero de elemento total. $\alpha_{add}^{th}$ es la sensibilidad del elemento clasificado justo debajo del ultimo elemento añadido. $\alpha_{add}^{th}$ es determinado para que el volumen eliminado sea igual a $V_k - V_{k-1} + volumen del elemento añadido$.

  
El ciclo de análisis de elemento finitos y remover/añadir elementos continua hasta que se llegue a la restricción de volumen $V^*$ o se consiga el siguiente criterio de convergencia

$$
\begin{equation}
error = \frac{\sum^N_{i=1} C_{k-i+1} - \sum^N_{i=1} C_{k-N-i+1}}{\sum^N_{i=1} C_{k-i+1}} \leq \tau
\end{equation}
$$

  

donde $k$ es la iteración actual, $\tau$ es la tolerancia de la convergencia permitida y $N$ es un numero entero. Normalmente, $N$ es $5$, ya que implica que el cambio en el cumplimiento mínimo sobre las ultimas 10 iteraciones es aceptablemente pequeño.

  
El método BESO puede ser resumido en los siguientes pasos:

1.  Discretizar el dominio y asignar la propiedad inicial, 0 o 1, para cada elemento.
2.  Realizar el análisis de elementos finito y calcular la sensibilidad de acuerdo con (13).
3.  Promediar la sensibilidad de la sensibilidad con el historial y guardar el resultado para la siguiente iteración (15).
4.  Determinar el volumen objetivo para la siguiente iteración (16).
5.  Añadir y eliminar elementos.
6.  repetir los pasos del 2 al 5 hasta llegar al volumen $V^*$ o que el criterio de convergencia (17) se satisfaga.

## Referencias
- Yago, D., Cante, J., Lloberas-Valls, O., & Oliver, J. (2022). Topology Optimization Methods for 3D Structural Problems: A Comparative Study. Archives of Computational Methods in Engineering, 29(3), 1525–1567. doi:10.1007/s11831-021-09626-2
- Huang, X., & Xie, M. (2010). Evolutionary Topology Optimization of Continuum Structures: Methods and Applications. John Wiley & Sons.

