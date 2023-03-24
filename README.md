# A Bayesian model of motion extrapolation explains offset mislocalisations with dynamic noise

To reproduce the figures from the paper, install the encessary packages with,

```
julia> using Pkg
julia> Pkg.activate("path/to/this/project")
julia> Pkg.instantiate()
```

and then run the file `main.jl`.


Estimation procedure            |  Estimated position
:-------------------------:|:-------------------------:
![A gif showing the particles' estimated position, overlayed to the illusion](https://github.com/phijoh/kalman-project/blob/master/readme-figures/localisation.gif)  |  ![A time series of the estimated position in the case of static and dynamic noise](https://github.com/phijoh/kalman-project/blob/master/readme-figures/position.png)
