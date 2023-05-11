---
title: "ESO/BESO"
date: 2023-03-15T11:48:52-05:00
draft: false
math: true
---

The ESO (evolutionary structural optimization) and BESO (bi-directional evolutionary structural optimization) methods are methods with hard-kill optimization criteria, that is, they can eliminate material and this causes the domain to be modified.

The ESO method consists of the simple concept of removing material that is inefficient. This is done by calculating the stresses in the entire model and eliminating the elements with low values.

## Optimization by the level of stress

For this, the stress level of each element is calculated and compared with criteria defined as rejection criteria.

  

$$
\begin{equation}
	\frac{\sigma_e}{\sigma_m} < RR_i 
\end{equation}
$$

  
$\sigma_e$ is the Von Mises stress of the element, $\sigma_m$ is the Von Mises stress of the model, and $RR_i$ is the rejection criterion. This inequality defines which elements are removed after each finite analysis. This cycle is repeated until a state of equilibrium is reached, that is, there are no more elements to be removed with the rejection criteria used. After this, the value of $RR_i$ can be increased to capture elements with a higher level of stress and thus remove more material.

$$
\begin{equation}
	RR_{i+1} = RR_i + ER
\end{equation}
$$


$ER$ is an evolutionary rate that is added to the rejection criteria. The ESO method can be summarized as:

1. Discretization of the model in finite elements.
2. Perform a finite element analysis of the structure.
3. Remove elements that satisfy (1).
4. Increase the rejection criterion using (2).
5. Repeat steps 3 and 4 until you reach the desired point.


## Strain energy optimization.

This criterion consists of taking stiffness as a factor to carry out the optimization. The operation is not allowed because the subscription has exceeded its free quota. Commonly called mean compliance $C$, it is the inverse measure of the stiffness of the entire structure. The minimum compliance can be defined as the total strain energy

$$
\begin{equation}
C = \frac{1}{2} f^T u
\end{equation}
$$

where $f$ and $u$ are the vectors of force and displacements respectively. The equilibrium of the structure is defined by:

$$
\begin{equation}
Ku = f \end{equation}
$$
where $K$ is the global stiffness matrix. When an element is removed, the stiffness matrix is rewritten as follows:

$$K^* = K - K_i$$ $$\Delta K = K^* - K = -K_i$$


$K^*$ is the array modified after removing an element, $K_i$ is the array of the element. In the equation, after eliminating an element, only the displacement vector varies, which is described

$$
\begin{equation}
Δu = -K^{-1} ΔKu
\end{equation}
$$

from equations (3) and (5)

  

$$
\begin{equation}
\Delta C = \frac{1}{2} f^T \Delta u = - \frac{1}{2} f^T k^{-1} \Delta K u = \frac{1}{2} u^T_i K_iu_i
\end{equation}
$$

  

where $u_i$ is the displacement vector of the element. $\Delta C$ is understood as the sensitivity number of minimum compliance and is denoted by $\alpha^e_i$.

  
$$
\begin{equation}
\alpha^e_i = \frac{1}{2}u^T_i K_iu_i
\end{equation}
$$


This equation expresses that the increase in strain energy or minimum compliance is equal to the strain energy of the element. To minimize the minimum compliance, which is equivalent to maximizing the stiffness matrix, elements with smaller values of $\alpha^e_i$ are eliminated.

ESO with stiffness optimization can be summarized in the following steps:

1. Discretization of the model in finite elements.
2. Perform a finite element analysis of the structure.
3. Calculate the sensitivity value of each element with (7).
4. Remove the elements with the lowest sensitivity value according to the rejection ratio.
5. Repeat steps 2 through 4 until minimum compliance or maximum strain is achieved.

The BESO method is an improved version of the ESO method, since this, in addition to eliminating material, also adds. Solid elements with low sensitivity value are removed and empty elements with high sensitivity value are added and become solid elements. The sensitivity value of the empty elements is calculated by means of extrapolation. For this method, the rejection ratio $RR$ defines the number of elements that are eliminated in each iteration and the inclusion ratio $IR$ is the one that defines the increase in the rejection ratio.

A problem that occurs in the ESO method is the high sensitivity it has when applying the restrictions, making it possible to achieve a worse design than the initial one. Therefore, in this new method, the optimization problem is better posed so that in this way better convergence is ensured.

The problem statement addressed by this method is defined by:

$$\begin{equation}
Minimize\ C = \frac{1}{2} f^T u
\end{equation}$$
$$\begin{equation}
V^* - \sum_{i=1}^N V_i x_i = 0
\end{equation}$$

  
$V^*$ is the volume of the entire structure and $v_i$ is the volume of the individual element. $x_i$ is a variable that defines when an element is present and when it is not, 1 or 0 respectively.

For this method, the sensitivity number is the one described by expression (7). When the mesh is non-uniform, the sensitivity number must take into account the volume of the element. In this way, the sensitivity number can be replaced by the strain energy density of the element.


$$
\begin{equation}
\alpha^e_i = \left(\frac{1}{2} u^T_i K_iu_i\right)/V_i
\end{equation}
$$

  
When a structure is discretized using a few elements, convergence problems can be generated because the sensitivity number becomes discontinuous. To avoid this problem, a filter was devised that softens the sensitivity by means of nodes and elements to avoid mesh dependency.

To carry out the smoothing filter, the nodes must first be smoothed, for this, the following equation is used
  

$$
\begin{equation}
\alpha^n_j = \sum_{i=1}^M w_i\alpha^e_i
\end{equation}
$$
  
where $M$ is the number of elements connected to node $j$. $w_i$ is the weight of element $i$ and $\sum_{i=1}^M wi = 1$.
  

$$
\begin{equation}
w_i = \frac{1}{M-1} \left( 1 - \frac{r_{ij}}{\sum_{i=1}^M r_{ij}} \right)
\end{equation}
$$

  

where $r_{ij}$ is the distance between the center of the element $i$ and the node $j$, in this way, the weight is greater when the distance between the element and the node is smaller. After performing the smoothing of the nodes, the smoothing of the elements is performed. To carry out this scheme, the distance $r_{min}$ is considered, which does not change with the refinement of the mesh. This distance defines the radius of the $\omega_i$ domain. The purpose of this domain is to identify the nodes that will influence the sensitivity of the element $i$.
  

$$
\begin{equation}
\alpha_i = \frac{\sum_{i=1}^K w(r_{ij}) \alpha_j^n}{\sum_{i=1}^K w(r_{ij})}
\end{equation}
$$

where K is the total number of nodes in the domain $\omega_i$, $w(r_{ij})$ is the defined weight factor

  

$$
\begin{equation}
w(r_{ij}) = r_{min} - r_{ij}\ (j=1,2,...,K)
\end{equation}
$$

In this way, the sensitivity number is smoothed and the sensitivity of the empty elements is obtained. With this method, it is possible to solve the problem generated by the dependency of the mesh, however, this generates large oscillations in the objective function, producing convergence problems. For this, it was found that averaging the current sensitivity with the previous one solves this problem, this average is described by means of
  
$$
\begin{equation}
\alpha_i = \frac{\alpha_i^k + \alpha_i^{k-1}}{2}
\end{equation}
$$

  

where $k$ is the current iteration.
To remove and add elements, the target volume of the next iteration $V_{k+1}$ must be previously calculated. Since the volume constraint $v^*$ can be higher or lower than the initial volume, the target volume can increase or decrease. The volume evolves as follows
  

$$
\begin{equation}
V_{k+1} = V_k(1 ± ER)\ (k=1,2,3,...)
\end{equation}
$$

  

where $ER$ is the volume evolution ratio. Once the volume constraint is met, the volume of the next iteration will remain constant.

When calculating the sensitivity of the elements, they are ordered from highest to lowest. For the solid elements that will be eliminated

  

$$\alpha_i \leq \alpha_{del}^{th}$$

  

for empty elements to be added
  

$$\alpha_i > \alpha_{add}^{th}$$

  

To calculate the values of $\alpha_{add}^{th}$ and $\alpha_{del}^{th}$, follow these steps:

  
1. Let $\alpha_{add}^{th} = \alpha_{del}^{th} = \alpha_{th}$. If there are 1000 elements and with the ordered sensitivities $\alpha_1 > \alpha_2 > ... > \alpha_{1000}$ and $V_{k+1}$ correspond to a design with 600 elements, then $\alpha_{th} = α_{600}$
  

2. Calculate the volume addition ratio $AR$, which is calculated by dividing the number of items added by the total number of items. If $AR \leq AR_{max}$ then skip step 3. Otherwise you have to recalculate $\alpha_{add}^{th}$ and $\alpha_{del}^{th}$.

  

3. To compute $\alpha_{add}^{th}$ sort the sensitivity of empty elements. The number of elements that will change from 0 to 1 will be equal to $AR_{max}$ multiplied by the number of total elements. $\alpha_{add}^{th}$ is the sensitivity of the element ranked just below the last element added. $\alpha_{add}^{th}$ is determined so that the volume removed is equal to $V_k - V_{k-1} + volume of the element added$.

  
The cycle of finite element analysis and remove/add elements continues until the volume constraint $V^*$ is reached or the next convergence criterion is met

$$
\begin{equation}
error = \frac{\sum^N_{i=1} C_{k-i+1} - \sum^N_{i=1} C_{k-N-i+1}}{\sum^N_{i=1 } C_{k-i+1}} \leq \tau
\end{equation}
$$

  

where $k$ is the current iteration, $\tau$ is the allowed convergence tolerance, and $N$ is an integer. Normally, $N$ is $5$, since it implies that the change in minimum compliance over the last 10 iterations is acceptably small.

  
The BESO method can be summarized in the following steps:

1. Discretize the domain and assign the initial property, 0 or 1, for each element.
2. Perform the finite element analysis and calculate the sensitivity according to (13).
3. Average the sensitivity of the sensitivity with history and save the result for the next iteration (15).
4. Determine the target volume for the next iteration (16).
5. Add and remove items.
6. Repeat steps 2 to 5 until the volume $V^*$ is reached or the convergence criterion (17) is satisfied.

## References
- Yago, D., Cante, J., Lloberas-Valls, O., & Oliver, J. (2022). Topology Optimization Methods for 3D Structural Problems: A Comparative Study. Archives of Computational Methods in Engineering, 29(3), 1525–1567. doi:10.1007/s11831-021-09626-2
- Huang, X., & Xie, M. (2010). Evolutionary Topology Optimization of Continuum Structures: Methods and Applications. John Wiley & Sons.