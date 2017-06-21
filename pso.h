/* An implementation of the Particle Swarm Optimization algorithm

   Copyright 2010 Kyriakos Kentzoglanakis

   This program is free software: you can redistribute it and/or
   modify it under the terms of the GNU General Public License version
   3 as published by the Free Software Foundation.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see
   <http://www.gnu.org/licenses/>.
*/


#ifndef PSO_H_
#define PSO_H_


// CONSTANTS
#define PSO_MAX_SIZE 100 // max swarm size
#define PSO_INERTIA 0.7298 // default value of w (see clerc02)


// === NEIGHBORHOOD SCHEMES ===

// global best topology
#define PSO_NHOOD_GLOBAL 0

// ring topology
#define PSO_NHOOD_RING 1

// Random neighborhood topology
// **see http://clerc.maurice.free.fr/pso/random_topology.pdf**
#define PSO_NHOOD_RANDOM 2



// === INERTIA WEIGHT UPDATE FUNCTIONS ===
#define PSO_W_CONST 0
#define PSO_W_LIN_DEC 1



// PSO SOLUTION -- Initialized by the user
typedef struct {

  double error;
  double *gbest; // should contain DIM elements!!

} pso_result_t;



// OBJECTIVE FUNCTION TYPE
typedef double (*pso_obj_fun_t)(double *, int, void *);



// PSO SETTINGS
typedef struct {

  int dim; // problem dimensionality
  double x_lo; // lower range limit
  double x_hi; // higher range limit
  double goal; // optimization goal (error threshold)

  int size; // swarm size (number of particles)
  int print_every; // ... N steps (set to 0 for no output)
  int steps; // maximum number of iterations
  int step; // current PSO step
  double c1; // cognitive coefficient
  double c2; // social coefficient
  double w_max; // max inertia weight value
  double w_min; // min inertia weight value

  int clamp_pos; // whether to keep particle position within defined bounds (TRUE)
  // or apply periodic boundary conditions (FALSE)
  int nhood_strategy; // neighborhood strategy (see PSO_NHOOD_*)
  int nhood_size; // neighborhood size
  int w_strategy; // inertia weight strategy (see PSO_W_*)

  void *rng; // pointer to random number generator (use NULL to create a new RNG)
  long seed; // seed for the generator

} pso_settings_t;



// return the swarm size based on dimensionality
int pso_calc_swarm_size(int dim);

// set the default PSO settings
void pso_set_default_settings(pso_settings_t *settings);


// minimize the provided obj_fun using PSO with the specified settings
// and store the result in *solution
void pso_solve(pso_obj_fun_t obj_fun, void *obj_fun_params,
	       pso_result_t *solution, pso_settings_t *settings);



#endif // PSO_H_
