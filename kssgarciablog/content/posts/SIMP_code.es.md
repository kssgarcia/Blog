---
title: "Código de optimización SIMP"
date: 2023-05-17T11:48:52-05:00
draft: false
math: true
---

Esta publicación presentará un código implementado en Python para el método de optimización topologica SIMP.

Primero, consideremos los siguientes [pasos](https://kssgarcia.github.io/posts/simp/) para hacer el algoritmo:

1. Discretizar el dominio e inicializar las variables de diseño.
2. Realice el análisis de elementos finitos.
3. Calcular la sensibilidad de la función objetivo con respecto a la densidad de cada elemento.
4. Realizar los criterios de optimalidad.
5. Calcular el error de la densidad.
6. Repita los pasos 2 a 5 hasta que se logre la masa $M^*$ o la convergencia del error.

``` python
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import animation
import numpy as np
from beams import *
from SIMP_utils import *

# Mesh
length = 50
height = 20
nx = 100
ny= 40
nodes, mats, els, loads, BC = beam(L=length, H=height, nx=nx, ny=ny, n=3)

# Calculate centers and volumes
niter = 60
centers = center_els(nodes, els)
Vi = volume(els, length, height, nx, ny)

# Initialize the design variables
r_min = np.linalg.norm(nodes[0,1:3] - nodes[1,1:3]) * 4
penal = 3
Emin=1e-9
Emax=1.0
volfrac = 0.5
change = 10
g = 0

# Initialize the density, sensitivity and the iteration history
rho = volfrac * np.ones(ny*nx,dtype=float)
sensi_rho = np.ones(ny*nx)
rho_old = rho.copy()
d_c = np.ones(ny*nx)
d_v = np.ones(ny*nx)
rho_data = []

for i in range(niter):

    if change < 0.01:
        print('Convergence reached')
        break

    mats[:,2] = Emin+rho**penal*(Emax-Emin)

    IBC, UG, rhs_vec = preprocessing(nodes, mats, els, loads)
    UC, *_ = postprocessing(nodes, mats[:,:2], els, IBC, UG, strain_sol=False)

    # Sensitivity analysis
    sensi_rho[:] = sensitivity_els(nodes, mats, els, UC, nx, ny)
    obj = ((Emin+rho**penal*(Emax-Emin))*sensi_rho).sum()
    d_c[:] = (-penal*rho**(penal-1)*(Emax-Emin))*sensi_rho
    d_v[:] = np.ones(ny*nx)
    d_c[:] = density_filter(centers, r_min, rho, d_c)

    # Optimality criteria
    rho_old[:] = rho
    rho[:], g = optimality_criteria(nx, ny, rho, d_c, d_v, g)

    # Compute the change
    change = np.linalg.norm(rho.reshape(nx*ny,1)-rho_old.reshape(nx*ny,1),np.inf)

```

## Funciones necesarias

``` python
def beam_2(L=10, H=10, F=-1000000, E=206.8e9, v=0.28, nx=20, ny=20):
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
    x, y, els = pre.rect_grid(L, H, nx, ny)
    mats = np.zeros((els.shape[0], 3))
    mats[:] = [E,v,1]
    nodes = np.zeros(((nx + 1)*(ny + 1), 5))
    nodes[:, 0] = range((nx + 1)*(ny + 1))
    nodes[:, 1] = x
    nodes[:, 2] = y
    mask_1 = (x == L/2) & (y > H/4)
    mask_2 = (x == L/2) & (y < -H/4)
    mask = np.bitwise_or(mask_1, mask_2)
    nodes[mask, 3:] = -1

    mask_loads = (x == -L/2) & (y < H/6) & (y > -H/6)
    loads_nodes = nodes[mask_loads, 0]
    loads = np.zeros((len(loads_nodes), 3))
    loads[:, 0] = loads_nodes
    loads[:, 2] = F
    BC = nodes[mask, 0]
    return nodes, mats, els, loads, BC

def sparse_assem(elements, mats, nodes, neq, assem_op, uel=None):
    """
    Assembles the global stiffness matrix
    using a sparse storing scheme

    The scheme used to assemble is COOrdinate list (COO), and
    it converted to Compressed Sparse Row (CSR) afterward
    for the solution phase [1]_.

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    assem_op : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    uel : callable function (optional)
      Python function that returns the local stiffness matrix.

    Returns
    -------
    kglob : sparse matrix (float)
      Array with the global stiffness matrix in a sparse
      Compressed Sparse Row (CSR) format.

    References
    ----------
    .. [1] Sparse matrix. (2017, March 8). In Wikipedia,
        The Free Encyclopedia.
        https://en.wikipedia.org/wiki/Sparse_matrix

    """
    rows = []
    cols = []
    stiff_vals = []
    mass_vals = []
    nels = elements.shape[0]
    for ele in range(nels):
        kloc, mloc = ass.retriever(elements, mats, nodes, ele, uel=uel)
        if mats.shape[1] == 3:
            kloc *= mats[elements[ele, 0], 2]
        ndof = kloc.shape[0]
        dme = assem_op[ele, :ndof]
        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        rows.append(glob_row)
                        cols.append(glob_col)
                        stiff_vals.append(kloc[row, col])
                        mass_vals.append(mloc[row, col])

    stiff = coo_matrix((stiff_vals, (rows, cols)), shape=(neq, neq)).tocsr()
    mass = coo_matrix((mass_vals, (rows, cols)), shape=(neq, neq)).tocsr()
    return stiff, mass

def is_equilibrium(nodes, mats, els, loads):
    """
    Check if the system is in equilibrium
    
    Get from: https://github.com/AppliedMechanics-EAFIT/SolidsPy/blob/master/solidspy/solids_GUI.py
    
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
    assem_op, bc_array, neq = ass.DME(nodes[:, -2:], els, ndof_el_max=8)
    stiff_mat, _ = sparse_assem(els, mats, nodes, neq, assem_op)

    rhs_vec = ass.loadasem(loads, bc_array, neq)
    disp = sol.static_sol(stiff_mat, rhs_vec)
    if not np.allclose(stiff_mat.dot(disp)/stiff_mat.max(), rhs_vec/stiff_mat.max()):
        equil = False

    return equil
    
def preprocessing(nodes, mats, els, loads):
    """
    Compute IBC matrix and the static solve.
    
    Get from: https://github.com/AppliedMechanics-EAFIT/SolidsPy/blob/master/solidspy/solids_GUI.py
    
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
        Static displacement solve
    rh_vec : ndarray 
        Vector of loads
    """   

    assem_op, bc_array, neq = ass.DME(nodes[:, -2:], els, ndof_el_max=8)

    # System assembly
    stiff_mat, _ = sparse_assem(els, mats, nodes[:, :3], neq, assem_op)
    rhs_vec = ass.loadasem(loads, bc_array, neq)

    # System solution
    disp = sol.static_sol(stiff_mat, rhs_vec)
    if not np.allclose(stiff_mat.dot(disp)/stiff_mat.max(),
                       rhs_vec/stiff_mat.max()):
        print("The system is not in equilibrium!")
    return bc_array, disp, rhs_vec


def postprocessing(nodes, mats, els, bc_array, disp, strain_sol=True):
    """
    Compute the nodes displacements, strains and stresses.
    
    Get from: https://github.com/AppliedMechanics-EAFIT/SolidsPy/blob/master/solidspy/solids_GUI.py
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    mats : ndarray
        Array with models materials
    els : ndarray
        Array with models elements
    bc_array : ndarray 
        Boundary conditions array
    disp : ndarray 
        Static solve
        
    Returns
    -------
    disp_complete : ndarray 
        Displacements at elements.
    strain_nodes : ndarray 
        Strains at elements
    stress_nodes : ndarray 
        Stresses at elements
    """   
    
    disp_complete = pos.complete_disp(bc_array, nodes, disp)
    strain_nodes, stress_nodes = None, None
    if strain_sol == True:  
        strain_nodes, stress_nodes = pos.strain_nodes(nodes, els, mats, disp_complete)
    
    return disp_complete, strain_nodes, stress_nodes


def optimality_criteria(nelx, nely, rho, d_c, d_v, g):
    """
    Optimality criteria method.

    Parameters
    ----------
    nelx : int
        Number of elements in x direction.
    nely : int
        Number of elements in y direction.
    rho : ndarray
        Array with the density of each element.
    d_c : ndarray
        Array with the derivative of the compliance.
    d_v : ndarray
        Array with the derivative of the volume.
    g : float
        Volume constraint.

    Returns
    -------
    rho_new : ndarray
        Array with the new density of each element.
    gt : float
        Volume constraint.

    """
    l1=0
    l2=1e9
    move=0.2
    rho_new=np.zeros(nelx*nely)
    while (l2-l1)/(l1+l2)>1e-3: 
        lmid=0.5*(l2+l1)
        rho_new[:]= np.maximum(0.0,np.maximum(rho-move,np.minimum(1.0,np.minimum(rho+move,rho*np.sqrt(-d_c/d_v/lmid)))))
        gt=g+np.sum((d_v*(rho_new-rho)))
        if gt>0 :
            l1=lmid
        else:
            l2=lmid
    return (rho_new,gt)


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
        Number of elements in x direction.
    ny : float
        Number of elements in y direction.

    Return 
    ----------
    V: float

    """

    dy = length / nx
    dx = height / ny
    V = dx * dy * np.ones(els.shape[0])

    return V

def sensitivity_els(nodes, mats, els, UC, nx, ny):
    """
    Calculate the sensitivity number for each element.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    mats : ndarray
        Array with models materials
    els : ndarray
        Array with models elements
    UC : ndarray
        Displacements at nodes
    nx : float
        Number of elements in x direction.
    ny : float
        Number of elements in y direction.

    Returns
    -------
    sensi_number : ndarray
        Sensitivity number for each element.
    """   
    sensi_number = np.zeros(els.shape[0])
    params = tuple(mats[1, :])
    elcoor = nodes[els[1, -4:], 1:3]
    kloc, _ = uel.elast_quad4(elcoor, params)

    sensi_number = (np.dot(UC[els[:,-4:]].reshape(nx*ny,8),kloc) * UC[els[:,-4:]].reshape(nx*ny,8) ).sum(1)

    return sensi_number

def density_filter(centers, r_min, rho, d_rho):
    """
    Performe the sensitivity filter.
    
    Parameters
    ----------
    centers : ndarray
        Array with the centers of each element.
    r_min : float
        Minimum radius of the filter.
    rho : ndarray
        Array with the density of each element.
    d_rho : ndarray
        Array with the derivative of the density of each element.
        
    Returns
    -------
    densi_els : ndarray
        Sensitivity of each element with filter
    """
    densi_els = np.zeros(centers.shape[0])
    for i, c in enumerate(centers):
        delta = centers-c
        r_ij = np.linalg.norm(delta, axis=1)
        H = r_min - r_ij
        H = np.maximum(0.0, H)
        densi_els[i] = (rho*H*d_rho).sum()/(H.sum()*np.maximum(0.001,rho[i]))
        
    return densi_els

def center_els(nodes, els):
    """
    Calculate the center of each element.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes.
    els : ndarray
        Array with models elements.
        
    Returns
    -------
    centers : 
        Centers of each element.
    """
    centers = np.zeros((els.shape[0], 2))
    for el in els:
        n = nodes[el[-4:], 1:3]
        center = np.array([n[1,0] + (n[0,0] - n[1,0])/2, n[2,1] + (n[0,1] - n[2,1])/2])
        centers[int(el[0])] = center

    return centers

def plot_mesh(elements, nodes, disp, E_nodes=None):
    """
    Plot contours for model

    Get from: https://github.com/AppliedMechanics-EAFIT/SolidsPy/blob/master/solidspy/solids_GUI.py

    Parameters
    ----------
    nodes : ndarray (float)
        Array with number and nodes coordinates:
         `number coordX coordY BCX BCY`
    elements : ndarray (int)
        Array with the node number for the nodes that correspond
        to each element.
    disp : ndarray (float)
        Array with the displacements.
    E_nodes : ndarray (float)
        Array with strain field in the nodes.

    """
    # Check for structural elements in the mesh
    struct_pos = 5 in elements[:, 1] or \
             6 in elements[:, 1] or \
             7 in elements[:, 1]
    if struct_pos:
        # Still not implemented visualization for structural elements
        print(disp)
    else:
        pos.plot_node_field(disp, nodes, elements, title=[r"$u_x$", r"$u_y$"],
                        figtitle=["Horizontal displacement",
                                  "Vertical displacement"])
        if E_nodes is not None:
            pos.plot_node_field(E_nodes, nodes, elements,
                            title=[r"",
                                   r"",
                                   r"",],
                            figtitle=["Strain epsilon-xx",
                                      "Strain epsilon-yy",
                                      "Strain gamma-xy"])

```

El código completo puede ser encontrado en el repositorio de [github](https://github.com/kssgarcia/OptTopolgy).


![Gif result](/SIMPani1.gif)
![Gif result](/SIMPani2.gif)

Estos son algunos de los resultados arrojados por el código para dos vigas sometidas a una carga vertical.


# Referencias

[1] Nabaki, K., Shen, J. y Huang, X. (2019). Optimización topológica evolutiva de estructuras continuas considerando falla por fatiga. _Materiales y diseño_, _166_, 107586. doi:10.1016/j.matdes.2019.107586