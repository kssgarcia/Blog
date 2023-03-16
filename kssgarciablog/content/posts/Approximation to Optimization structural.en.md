---
title: "Approach to structural optimization."
date: 2023-02-16T11:48:52-05:00
draft: false
---

This post will talk about structural optimization from a more general point of view, with the aim of strengthening the bases for future applications.

To understand in a simple way the background of structural optimization it is a good idea to separate the two terms and take the definition of each. Optimization consists of doing things better and a structure can be understood as an assembly of material that has the purpose of supporting loading. Therefore, structural optimization can be defined as the technique of making an assembly of material withstand loads in the best way [1]. This definition is an abstract idea of what structural optimization is, but it is a simple approach that assumes a more intuitive idea of the concept.  When it is mentioned to do the "best" of something usually what is sought is to minimize and maximize, a common case in engineering is that it seeks to reduce the expenses of a process, this translates to minimize the cost of operation.

Structural optimization uses a mathematical model, this consists of identifying an objective function, this function is the one that classifies the design and will be the result that will seek to minimize or maximize. Design variable, is a function or vector that defines the desire, can be geometry or material. State variable, represents the response of the structure, such as stresses, deformations, displacements, etc. 

Structural optimization can be divided into three cases:

- Size optimization: is the dimensions minimization and maximization.
- Shape optimization: refers to the domain e.i in a problem of values on the border shape optimization refers to choosing the domain that will make the structure present less effort.
- Topological optimization: it is the configuration of material, topological optimization allows the elimination of material creating holes in the model, this can be with the aim of minimizing the amount of material. This case is the most studied in engineering.

Displacements, stresses, deformation, weight, geometry, forces, rigidity are some of the variables that can be optimized. An optimization problem is to choose one of these variables as a function to minimize or maximize and use a different one as a constraint. 

### References

1. Peter W. Christensen, A. K. (2009). _An Introduction to Structural Optimization_. Retrieve of https://link.springer.com/book/10.1007/978-1-4020-8666-3
