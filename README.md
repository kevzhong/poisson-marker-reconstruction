# Poisson marker-function reconstruction

![Alt text](2d-example.png?raw=true "2d-example")
![Alt text](3d-example.png?raw=true "3d-example")


 An implementation of the Poisson solver for computing an indicator/marker-function or volume-of-fluid (VOF) field, $\phi_{ijk}$ from a Lagrangian front representation, characterised by normal vectors $\vec{n}$ and element size $\Delta s$.

 I have documented most of the details [here](https://kevzhong.github.io/poisson-marker-reconstruction/).

In CFD codes adopting front-tracking methods the main motivation for this computation arises when two-phase flows with variable fluid properties are considered. E.g. the density, viscosity, or thermal conductivity at an Eulerian cell $ijk$. Using density $\rho_{ijk}$ as an example, this would be computed as:

$$
\rho_{ijk} = \phi_{ijk}\rho_1 + (1-\phi_{ijk})\rho_2.
$$

In front-tracking methods however, the Eulerian field $\phi_{ijk}$ is not known, as we only have information about the Lagrangian front which we are tracking (typically a triangulated surface in 3D). The program in this repository fulfils this purpose by constructing the indicator field $\phi$, given a triangulated 3D Lagrangian mesh.

This methodology of constructing the indicator function from a Lagrangian front representation via a Poisson equation appears to have been mostly popularised by Tryggvason and co-workers (see references below) and is the formulation I have followed. In the computational geometry literature however, this present method appears as quite a well-known, popular algorithm known as **Poisson Surface Reconstruction** (Kazdahn et al., 2006) which is typically present in many meshing software packages, although there are some slight differences.

This implementation employs fast-Poisson solver methods (i.e. FFTs or DCTs), which presumes specific boundary conditions for the Eulerian domain. For the presently-considered single geometries, we generally have $\phi = 0$ uniformly at the edges of the bounding box, so using either of `DCT-II` or `FFT` is generally OK.  

Quick and dirty `MATLAB` implementation is found in `/matlab-src`. The Fortran implementation is contained in `/src`. Some starter/sample binary .stl files are contained in `/stl` for technical demonstration. It is assumed that the input .stl files will not contain any "problematic" features. I.e. the surface is orientable, with consistent normal vector orientation, does not contain feature edges and generally has a homogeneous triangle area distribution.

In the Fortran implementation, the indicator field is dumped as a binary file with an additional `.xmf` header file for convenient viewing in Paraview.


## References

 - G. Tryggvason, R. Scardovelli, S. Zaleski (2011). Direct numerical simulations of gas–liquid multiphase flows, Cambridge University Press.
 - D. Juric, G. Tryggvason (1998). Computations of boiling flows, Int. J. Multiph. Flow 24, 387–410.
 - G. Tryggvason, B. Bunner, A. Esmaeeli, D. Juric, N. Al-Rawahi, W. Tauber, J. Han, S. Nas, Y.-J. Jan (2001). A front-tracking method for the computations of multiphase flow, J. Comput Phys. 169, 708–759.
 - S. Shin, D. Juric (2002). Modeling three-dimensional multiphase flow using a level contour reconstruction method for front tracking without connectivity, J. Comput. Phys. 180, 427-470.
 - M. Kazhdan, M. Bolitho, H. Hoppe (2006). Poisson surface reconstruction, in Proc. 4th Eurographics Symp. Geom. Process, 7(4).