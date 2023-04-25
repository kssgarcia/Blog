---
title: "Código de optimización ESO"
date: 2023-04-23T11:48:52-05:00
draft: false
math: true
---

En esta publicación, presentaré un código Python para el método de optimización de ESO.

### Planteamiento del problema

Consideremos una viga en voladizo con una fuerza vertical de 10kN, dimensiones 20x10 metros, módulo joven de 206.8 GPa y relación de Poisson 0.28 figura 1 $RR=0.01$, $ER=0.005$.

![Scenario 1: Across columns](/ESO_model.png)
|:--:|
| <b>Figura 1. Diagrama de vigas.</b>|


Primero, consideremos los siguientes pasos para hacer el algoritmo:

1. Discretización del modelo en elementos finitos
2. Realizar el análisis de elementos finitos de la estructura.
3. Eliminar elementos con menos estrés de Von Mises.
4. Aumentar la tasa de rechazo.
5. Repita los pasos 2-4 hasta alcanzar la fracción de volumen.

### Discretización

La discretización se realiza con la siguiente función con la ayuda de numpy

``` python
def beam(L=10, H=10, F=-1000000, E=206.8e9, v=0.28, nx=20, ny=20):
    """
    Make the mesh for a cuadrilateral model.

    Parameters
    ----------
    L : float (optional)
        Beam's lenght
    H : float (optional)
        Beam's height
    F : float (optional)
        Vertical force.
    E : string (optional)
        Young module
    v : string (optional)
        Poisson ratio
    nx : int (optional)
        Number of element in x direction
    ny : int (optional)
        Number of element in y direction

    Returns
    -------
    nodes : ndarray
        Nodes array
    mats : ndarray (1, 2)
        Mats array
    els : ndarray
        Elements array
    loads : ndarray
        Loads array
    BC : ndarray
        Boundary conditions nodes

    """
    mats = np.array([[E,v]])
    x, y, els = pre.rect_grid(L, H, nx, ny)
    nodes = np.zeros(((nx + 1)*(ny + 1), 5))
    nodes[:, 0] = range((nx + 1)*(ny + 1))
    nodes[:, 1] = x
    nodes[:, 2] = y
    mask = (x==L/2)
    nodes[mask, 3:] = -1

    mask_loads = (x == -L/2) & (y < H/6) & (y > -H/6)
    loads_nodes = nodes[mask_loads, 0]
    loads = np.zeros((len(loads_nodes), 3))
    loads[:, 0] = loads_nodes
    loads[:, 2] = F
    BC = nodes[mask, 0]
    return nodes, mats, els, loads, BC
```

con esta función se definen las variables necesarias para hacer análisis de elementos finitos

``` python
length = 20
height = 10
nx = 50
ny= 20
nodes, mats, els, loads, BC = beam(L=length, H=height, nx=nx, ny=ny)
elsI,nodesI = np.copy(els), np.copy(nodes)
```

### Análisis de elementos finitos

Con estas variables ya podemos resolver el problema y obtener los desplazamientos, deformaciones y tensiones. Para calcular el análisis de elementos finitos se utiliza la biblioteca [Solidspy 1.0.16](https://github.com/AppliedMechanics-EAFIT/SolidsPy).

``` python 
def preprocessing(nodes, mats, els, loads):
    """
    Compute IBC matrix and the static solve.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    mats : ndarray
        Array with models materials
    els : ndarray
        Array with models elements
    loads : ndarray
        Array with models loads
        
    Returns
    -------
    bc_array : ndarray 
        Boundary conditions array
    disp : ndarray 
        Static displacement solve.
    """   

    assem_op, bc_array, neq = ass.DME(nodes[:, -2:], els, ndof_el_max=8)
    print("Number of elements: {}".format(els.shape[0]))

    # System assembly
    stiff_mat, _ = ass.assembler(els, mats, nodes[:, :3], neq, assem_op)
    rhs_vec = ass.loadasem(loads, bc_array, neq)

    # System solution
    disp = sol.static_sol(stiff_mat, rhs_vec)
    if not np.allclose(stiff_mat.dot(disp)/stiff_mat.max(),
                       rhs_vec/stiff_mat.max()):
        print("The system is not in equilibrium!")
    return bc_array, disp

def postprocessing(nodes, mats, els, IBC, UG):
    """
    Compute the nodes displacements, strains and stresses.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    mats : ndarray
        Array with models materials
    els : ndarray
        Array with models elements
    IBC : ndarray 
        Boundary conditions array
    UG : ndarray 
        Static solve.
        
    Returns
    -------
    disp_complete : ndarray 
        Displacements at elements.
    strain_nodes : ndarray 
        Strains at elements.
    stress_nodes : ndarray 
        Stresses at elements.
    """   

	disp_complete = pos.complete_disp(bc_array, nodes, disp)
	strain_nodes, stress_nodes = None, None
	strain_nodes, stress_nodes = pos.strain_nodes(nodes, els, mats, disp_complete)

	return disp_complete, strain_nodes, stress_nodes
```


### Optimización

Primero definamos el número de iteraciones `niter`, la relación de rechazo `RR`, la relación de evolución `ER` y la fracción de volumen óptima `V_opt`.

En la secuencia cíclica, primero verificamos si el sistema aún está en equilibrio; de lo contrario, el ciclo se detendrá; después de eso, resolvemos el problema y calculamos el estrés de Von Mises para cada elemento y lo escalamos. `mask_del` es un arreglo numpy de tipo bool con el mismo tamaño de `els`, las posiciones con valor `True` indicarán que el elemento en esa posición exacta será eliminado, con la función `protect_els()` evitaremos la remoción de los elementos donde se definieron las cargas y condiciones de contorno.

Después de borrar los elementos del arreglo tenemos que dar un paso más para asegurar que el análisis FEM no fallará con el nuevo arreglo `els`, para esto tenemos que restringir los nodos libres `del_node()`. Como último paso aumentamos la tasa de rechazo.

``` python


def is_equilibrium(nodes, mats, els, loads):
    """
    Check if the system is in equilibrium
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    mats : ndarray
        Array with models materials
    els : ndarray
        Array with models elements
    loads : ndarray
        Array with models loads
        
    Returns
    -------
    equil : bool
        Variable True when the system is in equilibrium and False when it doesn't
    """   

    equil = True
    DME, IBC, neq = ass.DME(nodes, els)
    KG = ass.assembler(els, mats, nodes, neq, DME)
    RHSG = ass.loadasem(loads, IBC, neq)
    UG = sol.static_sol(KG, RHSG)
    
    if not(np.allclose(KG.dot(UG)/KG.max(), RHSG/KG.max())):
        equil = False
    return equil

def strain_els(els, E_nodes, S_nodes):
    """
    Compute the elements strains and stresses.
    
    Parameters
    ----------
    els : ndarray
        Array with models elements
    E_nodes : ndarray
        Strains at nodes.
    S_nodes : ndarray
        Stresses at nodes.
        
    Returns
    -------
    E_els : ndarray (nnodes, 3)
        Strains at elements.
    S_els : ndarray (nnodes, 3)
        Stresses at elements.
    """   
    
    E_els = []
    S_els = []
    for el in els:
        strain_nodes = np.take(E_nodes, list(el[3:]), 0)
        stress_nodes = np.take(S_nodes, list(el[3:]), 0)
        strain_elemt = (strain_nodes[0] + strain_nodes[1] + strain_nodes[2] + strain_nodes[3]) / 4
        stress_elemt = (stress_nodes[0] + stress_nodes[1] + stress_nodes[2] + stress_nodes[3]) / 4
        E_els.append(strain_elemt)
        S_els.append(stress_elemt)
    E_els = np.array(E_els)
    S_els = np.array(S_els)
    
    return E_els, S_els

def volume(els, length, height, nx, ny):
    """
    Volume calculation.
    
    Parameters
    ----------
    els : ndarray
        Array with models elements.
    length : ndarray
        Length of the beam.
    height : ndarray
        Height of the beam.
    nx : float

    ny : float
    Return 
    ----------
    V: float

    """

    dy = length / nx
    dx = height / ny
    V   = dx * dy * els.shape[0]

    return V

def del_node(nodes, els):
    """
    Retricts nodes dof that aren't been used.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    els : ndarray
        Array with models elements

    Returns
    -------
    """   
    n_nodes = nodes.shape[0]
    for n in range(n_nodes):
        if n not in els[:, -4:]:
            nodes[n, -2:] = -1

def protect_els(els, loads, BC):
    """
    Compute an mask array with the elements that don't must be deleted.
    
    Parameters
    ----------
    els : ndarray
        Array with models elements
    loads : ndarray
        Array with models loads
    BC : ndarray 
        Boundary conditions nodes
        
    Returns
    -------
    mask_els : ndarray 
        Array with the elements that don't must be deleted.
    """   
    mask_els = np.ones_like(els[:,0], dtype=bool)
    protect_nodes = np.hstack((loads[:,0], BC)).astype(int)
    protect_index = None
    for p in protect_nodes:
        protect_index = np.argwhere(els[:, -4:] == p)[:,0]
        mask_els[protect_index] = False
        
    return mask_els

niter = 200
RR = 0.01
ER = 0.005
V_opt = volume(els, length, height, nx, ny) * 0.60


ELS = None
for _ in range(niter):

    if not is_equilibrium(nodes, mats, els, loads) or volume(els, length, height, nx, ny) < V_opt: break
    ELS = els
    
    IBC, UG = preprocessing(nodes, mats, els, loads)
    UC, E_nodes, S_nodes = postprocessing(nodes, mats, els, IBC, UG)
    E_els, S_els = strain_els(els, E_nodes, S_nodes)
    vons = np.sqrt(S_els[:,0]**2 - (S_els[:,0]*S_els[:,1]) + S_els[:,1]**2 + 3*S_els[:,2]**2)
    RR_el = vons/vons.max()
    mask_del = RR_el < RR
    mask_els = protect_els(els, loads, BC)
    mask_del *= mask_els
    
    els = np.delete(els, mask_del, 0)
    del_node(nodes, els)
    RR += ER
print(RR)
```

# Resultados

Como se puede observar en la _figura 2_ se presentan los resultados del ejercicio propuesto al inicio del post, para el gráfico _2.a_ se evidencia que el método sí converge a una estructura óptima, sin embargo se evidencian varios elementos libres que también pueden Se eliminan para obtener una mejor topología. Este problema se produce porque el método no tiene en cuenta la dependencia de la malla y el patrón de tablero de ajedrez. Para obtener mejores resultados, puede variar $RR$, $ER$ o también aumentar la cantidad de elementos finitos.

![Scenario 1: Across columns](/ESO_1.png)
|:--:|
| <b>Figura 2. Topologías ESO para vigas en voladizo con diferente fracción de volumen: (a) 80%; (b) 70%; c) 60 %; (d) 50%.</b>|

# Referencias

[1] Nabaki, K., Shen, J. y Huang, X. (2019). Optimización topológica evolutiva de estructuras continuas considerando falla por fatiga. _Materiales y diseño_, _166_, 107586. doi:10.1016/j.matdes.2019.107586

