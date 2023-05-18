---
title: "BESO optimization code"
date: 2023-04-24T11:48:52-05:00
draft: false
math: true
---

In this post, I will present a Python code for the BESO optimization method.

### Problem statement

Let us consider a cantilever beam with a vertical force of 10kN, dimensions 20x10 meters, a young module of 206.8 GPa and Poisson's ratio 0.28 figure 1 $RR=0.01$, $ERR=0.005$.

![Scenario 1: Across columns](/ESO_model.png)
|:--:|
| <b>Figure 1. Beams's diagram.</b>|


First, let's consider the following steps to do the algorithm:

1. Discretize the domain and assign the initial property, 0 or 1, for each element.
2. Perform the finite element analysis and calculate the sensitivity according to.
3. Average the sensitivity of the sensitivity with history and save the result for the next iteration.
4. Determine the target volume for the next iteration.
5. Add and remove items.
6. Repeat steps 2 to 5 until the optimal volume is reached or the convergence criterion is satisfied.

### Discretization

The discretization is done with the following function with the help of NumPy

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

with this function will define the necessary variables to do finite element analysis

``` python
length = 20
height = 10
nx = 50
ny= 20
nodes, mats, els, loads, BC = beam_2(L=length, H=height, nx=nx, ny=ny)
elsI,nodesI = np.copy(els), np.copy(nodes)

IBC, UG = preprocessing(nodes, mats, els, loads)
UCI, E_nodesI, S_nodesI = postprocessing(nodes, mats, els, IBC, UG)
```

### Finite element analysis

With these variables, we can now solve the problem and obtain the displacements, deformations and stresses. To compute the finite element analysis the library [Solidspy](https://github.com/AppliedMechanics-EAFIT/SolidsPy) is used.

``` python
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
    print("Number of elements: {}".format(els.shape[0]))

    # System assembly
    stiff_mat, _ = ass.assembler(els, mats, nodes[:, :3], neq, assem_op)
    rhs_vec = ass.loadasem(loads, bc_array, neq)

    # System solution
    disp = sol.static_sol(stiff_mat, rhs_vec)
    if not np.allclose(stiff_mat.dot(disp)/stiff_mat.max(),
                       rhs_vec/stiff_mat.max()):
        print("The system is not in equilibrium!")
    return bc_array, disp, rhs_vec


def postprocessing(nodes, mats, els, bc_array, disp):
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
    strain_nodes, stress_nodes = pos.strain_nodes(nodes, els, mats, disp_complete)
    
    return disp_complete, strain_nodes, stress_nodes
```


### Optimization
``` python
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
    stiff_mat, _ = ass.assembler(els, mats, nodes[:, :3], neq, assem_op)
    rhs_vec = ass.loadasem(loads, bc_array, neq)
    disp = sol.static_sol(stiff_mat, rhs_vec)
    if not np.allclose(stiff_mat.dot(disp)/stiff_mat.max(), rhs_vec/stiff_mat.max()):
        equil = False

    return equil

def protect_els(els, nels, loads, BC):
    """
    Compute an mask array with the elements that don't must be deleted.
    
    Parameters
    ----------
    els : ndarray
        Array with models elements
    nels : ndarray
        Number of elements
    loads : ndarray
        Array with models loads
    BC : ndarray 
        Boundary conditions nodes
        
    Returns
    -------
    mask_els : ndarray 
        Array with the elements that don't must be deleted.
    """   
    mask_els = np.zeros(nels, dtype=bool)
    protect_nodes = np.hstack((loads[:,0], BC)).astype(int)
    protect_index = None
    for p in protect_nodes:
        protect_index = np.argwhere(els[:, -4:] == p)[:,0]
        mask_els[els[protect_index,0]] = True
        
    return mask_els

def del_node(nodes, els, loads, BC):
    """
    Retricts nodes dof that aren't been used and free up the nodes that are in use.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    els : ndarray
        Array with models elements
    loads : ndarray
        Array with models loads
    BC : ndarray 
        Boundary conditions nodes

    Returns
    -------
    """   
    protect_nodes = np.hstack((loads[:,0], BC)).astype(int)
    for n in nodes[:,0]:
        if n not in els[:, -4:]:
            nodes[int(n), -2:] = -1
        elif n not in protect_nodes and n in els[:, -4:]:
            nodes[int(n), -2:] = 0


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

def sensitivity_els(nodes, mats, els, mask, UC):
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
    mask : ndarray
        Mask of optimal estructure
    UC : ndarray
        Displacements at nodes

    Returns
    -------
    sensi_number : ndarray
        Sensitivity number for each element.
    """   
    sensi_number = []
    for el in range(els.shape[0]):
        if mask[el] == False:
            sensi_number.append(0)
            continue
        params = tuple(mats[els[el, 2], :])
        elcoor = nodes[els[el, -4:], 1:3]
        kloc, _ = uel.elast_quad4(elcoor, params)

        node_el = els[el, -4:]
        U_el = UC[node_el]
        U_el = np.reshape(U_el, (8,1))
        a_i = 0.5 * U_el.T.dot(kloc.dot(U_el))[0,0]
        sensi_number.append(a_i)
    sensi_number = np.array(sensi_number)
    sensi_number = sensi_number/sensi_number.max()

    return sensi_number

def adjacency_nodes(nodes, els):
    """
    Create an adjacency matrix for the elements connected to each node.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes.
    els : ndarray
        Array with models elements.
        
    Returns
    -------
    adj_nodes : ndarray, nodes.shape[0]
        Adjacency elements for each node.
    """
    adj_nodes = []
    for n in nodes[:, 0]:
        adj_els = np.argwhere(els[:, -4:] == n)[:,0]
        adj_nodes.append(adj_els)
    return adj_nodes

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
    centers : ndarray, nodes.shape[0]
        Adjacency elements for each node.
    """
    centers = []
    for el in els:
        n = nodes[el[-4:], 1:3]
        center = np.array([n[1,0] + (n[0,0] - n[1,0])/2, n[2,1] + (n[0,1] - n[2,1])/2])
        centers.append(center)
    centers = np.array(centers)
    return centers

def sensitivity_nodes(nodes, adj_nodes, centers, sensi_els):
    """
    Calculate the sensitivity of each node.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    adj_nodes : ndarray
        Adjacency matrix of nodes
    centers : ndarray
        Array with center of elements
    sensi_els : ndarra
        Sensitivity of each element without filter
        
    Returns
    -------
    sensi_nodes : ndarray
        Sensitivity of each nodes
    """
    sensi_nodes = []
    for n in nodes:
        connected_els = adj_nodes[int(n[0])]
        if connected_els.shape[0] > 1:
            delta = centers[connected_els] - n[1:3]
            r_ij = np.linalg.norm(delta, axis=1) # We can remove this line and just use a constant because the distance is always the same
            w_i = 1/(connected_els.shape[0] - 1) * (1 - r_ij/r_ij.sum())
            sensi = (w_i * sensi_els[connected_els]).sum(axis=0)
        else:
            sensi = sensi_els[connected_els[0]]
        sensi_nodes.append(sensi)
    sensi_nodes = np.array(sensi_nodes)

    return sensi_nodes

def sensitivity_filter(nodes, centers, sensi_nodes, r_min):
    """
    Performe the sensitivity filter.
    
    Parameters
    ----------
    nodes : ndarray
        Array with models nodes
    sensi_nodes : ndarray
        Array with nodal sensitivity
    centers : ndarray
        Array with center of elements
    r_min : ndarra
        Minimum distance 
        
    Returns
    -------
    sensi_els : ndarray
        Sensitivity of each element with filter
    """
    sensi_els = []
    for i, c in enumerate(centers):
        delta = nodes[:,1:3]-c
        r_ij = np.linalg.norm(delta, axis=1)
        omega_i = (r_ij < r_min)
        w = 1/(omega_i.sum() - 1) * (1 - r_ij[omega_i]/r_ij[omega_i].sum())
        sensi_els.append((w*sensi_nodes[omega_i]).sum()/w.sum())
        
    sensi_els = np.array(sensi_els)
    sensi_els = sensi_els/sensi_els.max()

    return sensi_els

```

Firs, let's define the number of iteration `niter`, `ERR` the evolutionary rejection ratio and `t` the allowable convergence tolerance, the minimum radius `r_min` and the optimal volume `V_opt`. After defining these variables we compute the adjacency array for the nodes `adjacency_nodes()` and the center of each element `center_els()`, these arrays will be used to achieve the filter scheme.

The loop's sequence is define with:

``` python
niter = 200
ERR = 0.005
t = 0.01

r_min = np.linalg.norm(nodes[0,1:3] - nodes[1,1:3]) * 1.5
adj_nodes = adjacency_nodes(nodes, els)
centers = center_els(nodes, els)

Vi = volume(els, length, height, nx, ny)
V_opt = Vi.sum() * 0.50

# Initialize variables.
ELS = None
mask = np.ones(els.shape[0], dtype=bool)
sensi_I = None
C_h = np.zeros(niter)
error = 1000

for i in range(niter):
    # Calculate the optimal design array elements
    els_del = els[mask].copy()
    V = Vi[mask].sum()

    # Check equilibrium
    if not is_equilibrium(nodes, mats, els_del, loads): 
        print('Is not equilibrium')
        break
    
    # Storage the solution
    ELS = els_del

    # FEW analysis
    IBC, UG, rhs_vec = preprocessing(nodes, mats, els_del, loads)
    UC, E_nodes, S_nodes = postprocessing(nodes, mats, els_del, IBC, UG)

    # Sensitivity filter
    sensi_e = sensitivity_els(nodes, mats, els, mask, UC)
    sensi_nodes = sensitivity_nodes(nodes, adj_nodes, centers, sensi_e) #3.4
    sensi_number = sensitivity_filter(nodes, centers, sensi_nodes, r_min) #3.6

    # Avarage the sensitivity numbers to the historical information 
    if i > 0: 
        sensi_number = (sensi_number + sensi_I)/2 # 3.8
    sensi_number = sensi_number/sensi_number.max()

    # Check if the optimal volume is reached and calculate the next volume
    V_r = False
    if V <= V_opt:
        els_k = els_del.shape[0]
        V_r = True
        break
    else:
        V_k = V * (1 + ERR) if V < V_opt else V * (1 - ERR)

    # Remove/add threshold
    sensi_sort = np.sort(sensi_number)[::-1]
    els_k = els_del.shape[0]*V_k/V
    alpha_del = sensi_sort[int(els_k)]

    # Remove/add elements
    mask = sensi_number > alpha_del
    mask_els = protect_els(els[np.invert(mask)], els.shape[0], loads, BC)
    mask = np.bitwise_or(mask, mask_els)
    del_node(nodes, els[mask], loads, BC)

    # Calculate the strain energy and storage it 
    C = 0.5*rhs_vec.T@UG
    C_h[i] = C
    if i > 10: error = C_h[i-5:].sum() - C_h[i-10:-5].sum()/C_h[i-5:].sum()

    # Check for convergence
    if error <= t and V_r == True:
        print("convergence")
        break

    # Save the sensitvity number for the next iteration
    sensi_I = sensi_number.copy()
```

# Results

As can be seen in _figure 2_ the results of the exercise proposed at the beginning of the post are presented, in the graph _2.a_ it is evident that the method converges to an optimal structure. In this BESO program, the dependency of the mesh and the checkerboard pattern problems are considered under the filter scheme _2.d_.

![Scenario 1: Across columns](/BESO.png)
|:--:|
| <b>Figure 2. BESO topologies for cantilever bea with diferent volume's fraction: (a) 80%; (b) 70%; (c) 60%; (d) 50%.</b>|


# References

[1] Nabaki, K., Shen, J., & Huang, X. (2019). Evolutionary topological optimization of continuous structures considering fatigue failure. _Materials and design_, _166_, 107586. doi:10.1016/j.matdes.2019.107586