# Galaxy
The aim of this project is to construct (parts of) a hydrodynamic simulation
for galaxy evolution. 
The project didn't quite reach this, but hopefully the start here is complete enough.


# Theory/Background/Structure
See `background.pdf`

# Running

## Setup
On OSC, to run this code, first you need to make julia available (I tested on 1.8.5,
other version of julia may or may not throw interesting errors with this code, 
but I have not explored this)
```
module load julia/1.8.5
```
Now, to initialize the environment, start up the julia prompt by typing `julia`. Next, 
type `]` to enter package mode. The package prompt should look something like `(@v1.8) pkg >`, so finally type the two lines
```
activate .
instantiate
```

The project is now ready for use!
Typing `del` will leave the pkg prompt, and `exit()` or pressing `Ctrl-D` will leave the julia shell.

## Scripts

There are numerous scripts in the main directory of the project. Each script runs the simulation with different parameters/options to explore how well it is working so far. I know validation of a simulation would reqire more work, but unfortunantly I don't look into this here.
In general, scripts can be run with `julia scriptname.jl` or from the julia shell with `include("scriptname.jl")`. Running each script seperately is likely the safest option

- `two_body_test.jl` just simulates two point particles orbiting eachother. This is a test
of the numerical stability and verification gravity works as expected.
- `two_body_plot.jl` creates the associated plot
- `two_body_p.jl` is a variation of the two body code except where pressure/temperature are
the only physics and the particles are on a head-on-head collision.
- `static_eq.jl` tests if the simulation is stable with an isothermic sphere.
- `sedov.jl` tests the response of an explosion or shock wave from the origin.


# Dependencies
This package needs the following dependencies, however, these can all be installed automatically
with 
activate . + instantiate


- CSV
- DataFrames
- ForwardDiff
- Glob
- LinearAlgebra
- Logging
- NearestNeighbors
- Plots
- Printf
- QuadGK
- Random
- Roots
- StaticArrays
- StatsPlots
- TOML



