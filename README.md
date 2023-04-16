# Galaxy
The aim of this project is to construct a simple hydrodynamical simulation 
in julia.

The equations follow from ...


# Theory

I chose to use smoothed particle hydrodynamics (SPH) which is based on the lagrangian
particle formulation. I calculate forces in the frame of the particle, using lagrangian
derivatives, i.e.
$$
\frac{d/dt} \def \frac{\partial}{\partial t} + \vec{v} \cdot \nabla
$$
