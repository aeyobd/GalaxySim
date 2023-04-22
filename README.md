# Galaxy
The aim of this project is to construct a simple hydrodynamical simulation in julia. 


I chose julia because this language is both high-level and compiled, so writing code is much easier than in C/C++, but it doesn't have the same performance penalty as python.

On OSC, to run this code, first you need to load julia
```
module load julia/1.8.5
```


# Theory

I chose to use smoothed particle hydrodynamics (SPH) which is based on the lagrangian
particle formulation. I calculate forces in the frame of the particle, using lagrangian
derivatives, i.e.

$$\alpha + \int_0^1~$$

***

# Implementation

## Density Calculation

Following 

I use Newton's method to find the density. The ideal smoothing length is given by 
$$
h=\eta\left(\frac{m}{\rho}\right)^{1/3},
$$
however, this equation is dependent on density above. So, we can find self-consistant values of $h$ and $\rho$ by finding the zero of 
$$
f(h) = \rho\vert_{h} - \rho_{\rm new}
$$
Where the first term and second term compare the current value of $\rho$ assuming $h$, and the calculated $\rho$ by taking the sum above.


![](density.png)
An example solution of the equation.

