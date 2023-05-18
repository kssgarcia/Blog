---
title: "SIMP Method for Topology Optimization"
date: 2023-03-20T11:48:52-05:00
draft: false
math: true
---

This method has been shown to be used in a large number of problems, since it is computationally efficient and simple to implement. Unlike BESO [ESO/BESO](https://kssgarcia.github.io/es/posts/eso_beso/) this method uses a *soft-kill* criterion, that is, the elements are not completely eliminated, they are only reduces the elastic modulus or one of the dimensional characteristics. This criterion is implemented since the complete removal of elements can result in theoretical difficulties.

The problem presented seeks to minimize the minimum compliance of the model. This is given by

  
$$
\begin{equation}
\text{Minimize } C = f^T u
\tag{1a}
\end{equation}
$$

$$
\begin{equation}
\text{Subject to: } V^* - \sum _{i=1}^N v_i x_i = 0
\tag{1b}
\end{equation}
$$

$$
\begin{equation}
x_i = x_{min} \text{ o } 1
\tag{1c}
\end{equation}
$$


$v_i$ is the volume of the element, $V^*$ is the volume chosen for the final structure and $x_i$ is the design variable that denotes the relative density of the element $i$; This method uses twice the strain energy that is used in the [ESO/BESO](https://kssgarcia.github.io/en/posts/eso_beso/) method. Unlike the KISS method, the layout variable is never zero, so elements are never completely removed, avoiding having empty elements; for the elements with $x_i=x_{min}$ they are called soft elements or *soft-elements* in English, the value of the variable $x_{min}$ can be for example $0.001$ but never zero.

In this method a material interpolation scheme with penalty is used to direct the solution to almost solid designs. To achieve this type of design, the modulus of elasticity is interpolated as a function of the element density or the design variable:


$$
\begin{equation}
E(x_i) = E_1x_i^p
\tag{2}
\end{equation}
$$

  
where $E_1$ is the young's modulus of the solid material element and $p$ is the penalty exponent. It is assumed that the Poisson ratio is independent of the design variable and the global stiffness matrix is expressed by the element matrix and the design variable $x_i$ as
  

$$
\begin{equation}
K = \sum_i x_i^p K_i^0
\tag{3}
\end{equation}
$$

  
where $K^0_i$ is the stiffness matrix of the solid element.

### Sensitivity analysis

To obtain the sensitivity number, first consider the equilibrium equation that is expressed by


$$
\begin{equation}
Ku = f
\tag{4}
\end{equation}
$$


if the minimum compliance (1a) is considered and it is assumed that the design variable $x_i$ varies continuously between $x_{min}$ and $1$, the sensitivity of the objective function with respect to the design variable is
  

$$
\begin{equation}
\frac{dC}{dx_i} = \frac{df^T}{dx_i}u + f^T\frac{du}{dx_i}
\tag{5}
\end{equation}
$$
  

The adjoint method is used to calculate the sensitivity of the displacement vector. By introducing a vector of Lagrange multipliers $\lambda$, an extra term $\lambda^T(f-Ku)$ can be added to the objective function without changing anything due to equation (4)
  
$$
\begin{equation}
C = f^Tu + \lambda^T(f-Ku)
\tag{6}
\end{equation}
$$

taking into account the variation of the design variable
  

$$
\begin{equation}
\frac{dC}{dx_i} = \frac{df^T}{dx_i}u + f^T\frac{du}{dx_i} + \frac{\lambda^T}{dx_i}(f-Ku) + \lambda^T\left(\frac{f^T}{dx_i}-\frac{K}{dx_i}u - K\frac{u}{dx_i}\right)
\tag{7}
\end{equation}
$$


In the equation shown above, terms number 1 and 3 can be eliminated because it is assumed that the variation of $x_i$ has no effect on the force vector, therefore $\frac{df}{dx_i}$ and therefore that $f-Ku = 0$ due to (4)
  
$$
\begin{equation}
\frac{dC}{dx_i} = (f^T - \lambda^T K) \frac{du}{dx_i} - \lambda^T \frac{dK}{dx_i} u
\tag{8}
\end{equation}
$$


from equation (6) since $(f-Ku)$ is equal to zero, the vector $\lambda$ can be chosen as is. To eliminate the term $\frac{du}{dx_i}$ from the previous equation, $\lambda$ is defined taking into account


$$
\begin{equation}
f^T - \lambda^T K = 0
\tag{9}
\end{equation}
$$


With the previous equation and the equilibrium equation (4) the following solution is obtained for $\lambda$


$$
\begin{equation}
\lambda = u
\tag{10}
\end{equation}
$$


substituting $\lambda$ into equation (8),


$$
\begin{equation}
\frac{dC}{dx_i} = - u^T \frac{dK}{dx_i} u
\tag{11}
\end{equation}
$$


Substituting equation (3) in the previous equation, we get

  
$$
\begin{equation}
\frac{dC}{dx_i} = - p x_i^{p-1} u^T_i K_i^0 u_i
\tag{12}
\end{equation}
$$

  
The previous optimization problem can be solved by different approaches such as the optimality criterion. A scheme for updating the layout variable can be


$$
x_i^{K+1}= \begin{cases}\max \left(x_{\min }, x_i^K-m\right) & \text { if } x_i^K B_i^\eta \leq \max \left( x_{\min }, x_i^K-m\right) \\\ \min \left(1, x_i^K+m\right) & \text { if } \min \left(1, x_i^K+m\right) \leq x_i^K B_i^\eta \\\ x_i^K B_i^\eta & \text { otherwise }
\end{cases}
\tag{13}
$$

where $x_i^K$ denotes the design value at iteration $K$, $m$ is the positive limit motion, $\eta$ is the damping ratio (usually $0.5$) and $B_i$ is obtained from the optimality condition


$$
\begin{equation}
B_i = \lambda^{-1} p x_i^{p-1} u_i^T K_i^0 u_i
\tag{14}
\end{equation}
$$


where $\lambda$ is the Lagrange multiplier that can be determined using the bisection or Newton method. To ensure that the optimal design is mesh-independent and checkerboard-free, the following filter scheme is introduced that modifies the sensitivity of the elements


$$
\begin{equation}
\frac{\partial c}{\partial x_i}=\frac{1}{x_i \sum_{j=1}^N H_{i j}} \sum_{j=1}^N H_{i j} x_j \frac {\partial c}{\partial x_j}
\tag{15}
\end{equation}
$$

  
where $N$ is the total number of elements in the mesh and $H_{ij}$ is the weight factor given by
  

$$
\begin{equation}
H_{i j}=r_{\min }-r_{i j}, \quad \lbrace i \in N | r_{ij} \leq r_{min} \rbrace
\tag{16}
\end{equation}
$$

where $r_{ij}$ is the distance between the centers of the elements $i$ and $j$. The weight $H_{ij}$ is zero outside the circle of radius $r_{min}$.

1. Discretized the domain and initialize the design variables.
2. Perform the finite element analysis.
3. Calculate the sensitivity of the objective function with respect to the density of each element.
4. Perform the optimality criteria.
5. Compute the error of the density.
6. Repeat steps 2 to 5 until the mass $M^*$ or the convergence of the error is achieved. 

## References

- Yago, D., Cante, J., Lloberas-Valls, O., & Oliver, J. (2022). Topology Optimization Methods for 3D Structural Problems: A Comparative Study. Archives of Computational Methods in Engineering, 29(3), 1525–1567. doi:10.1007/s11831-021-09626-2

- Huang, X., & Xie, M. (2010). Evolutionary Topology Optimization of Continuum Structures: Methods and Applications. John Wiley & Sons.