---
mathjax: true
---
It is a FEM based realization of the method of boundary algebraic equations for Helmholtz equation on a square lattice.

![Sample geometry](some_geometry.png)

We consider a problem of diffraction by a finite scatterer with Dirichlet boundary condition imposed on the surface:

$$ u(m,n) = 0, \quad (m,n) \in \Omega_s $$,

where $\Omega_s$ is the set of surface nodes. 

The homogeneous discrete Helmholtz equation

$$ \Delta_{(m,n)}[u] + k^2 u_{(m,n)} = 0  $$

is satisfied in the rest of the domain. The operator $\Delta_{(m,n)}[\cdot]$ is the 2D discrete Laplace operator, which is a 5-point finite difference approximation of the continuos Laplace operator:

$$ \Delta_{(m,n)}[u] = u(m,n+1) + u(m,n-1) + u(m+1,n) + u(m-1,n)-4u(m,n). $$

By changing element matrices $K_{el}$ and $M_{el}$ in the script Assemble_K_M one can introduce more complicated approximations of the Laplace operator, as the 9-point FEM stencil.


