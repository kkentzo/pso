Particle Swarm Optimization (PSO) in C
===



An implementation of the Particle Swarm Optimization (PSO) algorithm
[1,2] in C that can be "plugged into" your code as a small
library. PSO is used for problems involving global stochastic
optimization of a continuous function (called the objective
function). PSO can also be used for discrete optimization problems,
but this behaviour is not implemented in the current version of this
library.

There's also an implementation in Go which can be found
[here](https://github.com/kkentzo/pso.go)







## USAGE


Just include pso.h and pso.c in your code. You need to have the [GNU
Scientific Library](http://www.gnu.org/software/gsl/) and the
respective development (i.e. header) files in order to include pso.c
in your application. In your Makefile add `-lgsl` and `-lgslcblas` to
your `LDFLAGS`.

In order to use `pso_solve()`, you need :

1. an objective function to be minimized (see defined type
`pso_obj_fun_t` in pso.h),

2. a `pso_results_t` object with a properly initialized (malloc'd)
gbest buffer. This is where the best position discovered will be
stored, along with the minimum error (stored in member `error`).

3. a `pso_settings_t` object with properly initialized values (use
`pso_set_default_settings()` for a quick and dirty initialization)






## FEATURES




### NEIGHBORHOOD TOPOLOGIES

pso provides three different strategies for determining each
particle's neighborhood attractor:

1. global topology (`PSO_NHOOD_GLOBAL`) where each particle is informed
by every other particle in the swarm

2. ring topology (`PSO_NHOOD_RING`) where there exists a fixed ring
topology and each particle is informed by its adjacent particles

3. random topology (`PSO_NHOOD_RANDOM`) where a random topology is
used that is updated when the error is not improved in two successive
iterations [3,4]. Use `nhood_size` in `pso_settings_t` to adjust the
average number of informers for each particle.




### INERTIA WEIGHT STRATEGIES

The value of the inertia weight (w) determines the balance between
global and local search. Two different strategies are implemented:

1. Constant value (PSO_W_CONST) using w=0.7298 [5]

2. Linearly decreasing inertia weight (`PSO_W_LIN_DEC`) [6]. Use
`w_max` and `w_min` in `pso_settings_t` to control the starting and
ending point respectively.




## OTHER SETTINGS


The following can be set in the `pso_settings_t` struct:


`dim` : problem dimensionality which should be equal to the size of
the gbest buffer in `pso_result_t` as well as to the size of the first
argument (pointer to a position buffer) of the objective function
(`pso_obj_fun_t`)

`size` : the number of particles in the swarm (the function
`pso_calc_swarm_size()` is also provided for the automatic calculation
of th swarm size based on the problem dimensionality).

`rng` : a pointer to a `gsl_rng` object (from GNU gsl). Pass NULL for
automatic generation of the random number generator.

`seed` : the seed to use for `rng` (default is time(0))

`x_lo, x_hi` : boundaries for particle positions

`clamp_pos` : if TRUE then the position of a particle that has exceeded
a boundary is set to the value of that boundary and its velocity
becomes 0 (along the dimension where the boundary was exceeded). If
it's FALSE, then periodic boundary conditions are enforced.

`print_every` : if greater than zero then this value specifies how many
steps should pass before information about the state of the search is
printed on screen

`c1, c2` : cognitive and social coefficients respectively. The default
values are c1 = c2 = 1.496 [5].

`steps` : the maximum number of steps to run the algorithm for.

`goal` : if the objective function returns a value lower than this
goal the search will stop








## EXAMPLES

A file demo.c with its Makefile are provided for your
convenience. demo.c provides instructions on how to setup pso in your
application.






## DISCLAIMER

Feel free to use the code as you see fit; I accept no responsibility
for any damage/catastrophy that might occur as a result :)







## REFERENCES


[1] Kennedy J and Eberhart R, "Particle Swarm Optimization."
Proc. IEEE Int’l Conf. Neural Networks, vol. 4, pp. 1942-1948, 1995.

[2] Poli R, Kennedy J and Blackwell T. "Particle swarm
optimization." Swarm intelligence 1.1 (2007): 33-57.

[3] Xu, R., Wunsch II, D., & Frank, R. (2007). Inference of genetic
regulatory networks with recurrent neural network models using
particle swarm optimization. IEEE/ACM Transactions on Computational
Biology and Bioinformatics (TCBB), 4(4), 681-692.

[4] http://clerc.maurice.free.fr/pso/random_topology.pdf

[5] Clerc, M., & Kennedy, J. (2002). The particle swarm-explosion,
stability, and convergence in a multidimensional complex
space. Evolutionary Computation, IEEE Transactions on, 6(1), 58-73.

[6] Shi, Y., & Eberhart, R. (1998). Parameter selection in particle
swarm optimization. In Evolutionary Programming VII
(pp. 591-600). Springer Berlin/Heidelberg.
