# Governing formulation

I follow mostly the working of Tryggvason et al. (2011, $\S$ 6.5.1).

The problem is as follows. Given an interface front located at $\mathbf{x} = \mathbf{x}_f$, we seek to determine the marker-function field $\phi(x,y,z)$ given discrete knowledge of the interface front. First, we can relate the jump in the indicator function (i.e. its gradient), $\mathbf{G}$ analytically via

$$
\mathbf{G} \equiv \nabla \phi = \int\Delta\phi \mathbf{\hat{n}}\delta(\mathbf{x} - \mathbf{x}_f)\,\mathrm{d}A
$$

where $\Delta \phi$ is the jump in the indicator function across the interface. Typical convention is to let $\phi \in [0,1]$ such that $\Delta \phi = 1$ as a constant (although this also depends on the convention adopted for the normal vector direction).

The Poisson reconstruction method seeks to construct the indicator field by solving a Poisson equation for $\phi$. That is, a divergence form of the above equation:

$$
\nabla^2 \phi = \nabla \cdot \mathbf{G}.
$$

# Numerical discretisation

Suppose the Lagrangian front $\mathbf{x}_f$ is represented by $l$ Lagrangian nodes, characterised by element areas $\Delta s_l$ and normal vectors $\mathbf{\hat{n}}_l$. 

To discretise the RHS of the Poisson equation, we employ a standard Lagrangian-to-Eulerian interpolation procedure:

$$
\mathbf{G}_{ijk} = \Delta \phi \sum_{l} \mathbf{\hat{n}}_l \tilde{\delta}^l_{ijk}\frac{\Delta s_l}{\Delta x \Delta y \Delta z}
$$

where $\tilde{\delta}^l_{ijk}$ is a regularised delta function. Again, common in Lagrangian methods. It is computed as a product of 1D regularised delta functions, $d(\cdot)$:

$$
\tilde{\delta}^l_{ijk} = d(r_x)d(r_y)d(r_z), \quad r_x = (x_l - x_i) / \Delta x, \, \text{etc}.
$$

Various choices exist for the regularised delta functions (see for example table 6.1 of Tryggvason et al. 2011, table 5.1 of Kajishima & Taira 2016). A choice which has worked well from my experience in testing is the formulation of Brackbill & Ruppel (1986):

$$
d(r) = \begin{cases}
\frac{2}{3} - |r|^2 + \frac{1}{2}|r|^3, \quad 0 \leq |r| \leq 1, \\
\frac{1}{6}(2- |r|)^3, \quad 1 \leq |r| \leq 2,\\
0, \qquad \text{otherwise}.
\end{cases}
$$

Next we consider the computation of $\mathbf{G}_{ijk}$. Notice how in the RHS of the Poisson equation, we want $\nabla \cdot \mathbf{G}$. Since our indicator field will be stored at cell centres, we ideally want $\nabla \cdot \mathbf{G}$ situated at cell faces. That is, we should employ a staggered grid calculation when evaluating $\mathbf{G}$ from the interpolation. Writing $\mathbf{G} = (G^x,G^y,G^z)$, then we should have for example...

$$
G^x_{i\pm1/2,j,k} = \Delta \phi \sum_{l} \hat{n}^x_l \tilde{\delta}^l_{i\pm1/2,j,k}\frac{\Delta s_l}{\Delta x \Delta y \Delta z}.
$$

With this, we can then write out the stencil of the discrete Poisson equation explicitly, $\nabla^2 \phi_{ijk} = (\nabla \cdot G)_{ijk}$:

$$
\frac{\phi_{i+1,j,k} - 2\phi_{i,j,k} + \phi_{i-1,j,k}}{\Delta x^2} + \frac{\phi_{i,j+1,k} - 2\phi_{i,j,k} + \phi_{i,j-1,k}}{\Delta y^2} + \frac{\phi_{i,j,k+1} - 2\phi_{i,j,k} + \phi_{i,j,k-1}}{\Delta z^2} = \frac{1}{\Delta x} \left[ G^x_{i+1/2,j,k} - G^x_{i-1/2,j,k} \right] + \frac{1}{\Delta y} \left[ G^y_{i,j+1/2,k} - G^x_{i,j-1/2,k} \right] + \frac{1}{\Delta z} \left[ G^z_{i,j,k+1/2} - G^x_{i,j,k-1/2} \right]
$$

This discrete Poisson equation can be solved using any standard method, but some notes on boundary conditions are in order. We could start the solution (applying Dirichlet conditions) from locations where we know the marker function (i.e. external points). To treat walls, we can also exclude the stencil of the wall-normal direction at cells adjacent to the wall.

In the present implementation, I have presumed the bounding Eulerian box is "uniform" such that I can apply uniform Dirichlet $\phi = 0$ for all boundaries. Periodic boundary conditions also work for this. This enables the use of fast-Poisson solver methods via DCTs or FFTs. That is, I solve:

$$
\left( \frac{\lambda_{k}}{\Delta x^2} + \frac{\lambda_{m}}{\Delta y^2} + \frac{\lambda_{n}}{\Delta z^2} \right)\hat{\phi}_{kmn} = \widehat{(\nabla \cdot \mathbf{G})}_{kmn} \implies \hat{\phi}_{kmn} = \frac{\widehat{(\nabla \cdot \mathbf{G})}_{kmn}}{\frac{\lambda_{k}}{\Delta x^2} + \frac{\lambda_{m}}{\Delta y^2} + \frac{\lambda_{n}}{\Delta z^2}}
$$

where $\lambda$ are the modified wavenumbers. This can be directly solved algebraically for each mode $(k,m,n)$ after which $\phi_{ijk}$ may be obtained via an inverse transform.

Some delicate consideration is needed when considering the $\lambda_k = \lambda_m = \lambda_n = 0$ mode to avoid a division by zero. Here, $\hat{\phi}_{000}$ corresponds to the mean value of $\phi_{ijk}$ in the entire domain, which physically, can be related to the volume of the immersed object, $\mathcal{V}$. Since we know the volume of the triangulated geometry (calculting the volume of a triangulated surface is a standard computation). If we have $\phi = 0$ outside the triangulated front and $\phi = 1$ inside, then we can set

$$
\hat{\phi}_{0 0 0}  =\mathcal{V} / (L_x L_y L_z)
$$

 to enforce volume conservation between the Eulerian indicator field and Lagrangian geometry. Here, $L_x L_y L_z$ is the total volume of the Eulerian domain. This is not without issues however, as this does not guarantee the boundedness of $\phi \in [0,1]$. From my experience, I have observed $\phi \in [\phi_{min},\phi_{max}] [-10^{-4},1.0008]$, close but not quite. This unboundedness can be catostraphic for multiphase CFD calculations, for instance by introducing negative densities, which are unphysical. A compromise is to shift the field from its bounds $[\phi_{min},\phi_{max}]$ to $\phi \in [0,1]$ using a uniform shift:

$$
\phi_{ijk} \leftarrow \frac{\phi_{ijk} - \phi_{min}}{\phi_{max} - \phi_{min}}.
$$

The trade-off, then, is a loss of volume conservation to ensure $\phi \in [0,1]$ boundedness. As of yet, I am unsure of any simple way in which both can be simultaneously satisfied.