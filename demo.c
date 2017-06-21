/* An implementation of the Particle Swarm Optimization algorithm

   === DEMONSTRATION CODE ===

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


#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include "pso.h"





//==============================================================
//                  BENCHMARK FUNCTIONS
//==============================================================

double pso_sphere(double *vec, int dim, void *params) {

  double sum = 0;
  int i;
  for (i=0; i<dim; i++)
    sum += pow(vec[i], 2);

  return sum;
}



double pso_rosenbrock(double *vec, int dim, void *params) {

  double sum = 0;
  int i;
  for (i=0; i<dim-1; i++)
    sum += 100 * pow((vec[i+1] - pow(vec[i], 2)), 2) +	\
      pow((1 - vec[i]), 2);

  return sum;

}


double pso_griewank(double *vec, int dim, void *params) {

  double sum = 0.;
  double prod = 1.;
  int i;
  for (i=0; i<dim;i++) {
    sum += pow(vec[i], 2);
    prod *= cos(vec[i] / sqrt(i+1));
  }

  return sum / 4000 - prod + 1;

}



//==============================================================

// BENCHMARK FUNCTION SETTINGS
void pso_set_sphere_settings(pso_settings_t *settings) {

  settings->x_lo = -100;
  settings->x_hi = 100;
  settings->goal = 1e-5;

}



void pso_set_rosenbrock_settings(pso_settings_t *settings) {

  settings->x_lo = -2.048;
  settings->x_hi = 2.048;
  settings->goal = 1e-5;

}

void pso_set_griewank_settings(pso_settings_t *settings) {

  settings->x_lo = -600;
  settings->x_hi = 600;
  settings->goal = 1e-5;

}



//==============================================================



int main(int argc, char **argv) {

  // define objective function
  pso_obj_fun_t obj_fun = pso_sphere;
  // initialize pso settings
  pso_settings_t settings;
  // set the default settings
  pso_set_default_settings(&settings);
  // set the problem specific settings
  pso_set_sphere_settings(&settings);

  // set PSO settings manually
  settings.size = 30;
  settings.nhood_strategy = PSO_NHOOD_RING;
  settings.nhood_size = 10;
  settings.w_strategy = PSO_W_LIN_DEC;

  // parse command line argument (function name)
  if (argc == 2) {
    if (strcmp(argv[1], "rosenbrock") == 0) {
      obj_fun = pso_rosenbrock;
      pso_set_rosenbrock_settings(&settings);
      printf("Optimizing function: rosenbrock (dim=%d, swarm size=%d)\n",
             settings.dim, settings.size);
    } else if (strcmp(argv[1], "griewank") == 0) {
      obj_fun = pso_griewank;
      pso_set_griewank_settings(&settings);
      printf("Optimizing function: griewank (dim=%d, swarm size=%d)\n",
             settings.dim, settings.size);
    } else if (strcmp(argv[1], "sphere") == 0) {
      printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
             settings.dim, settings.size);
    }
  } else if (argc > 2) {
    printf("Usage: demo [PROBLEM], where problem is optional with values [sphere|rosenbrock|griewank]\n ");
    return 1;
  }


  // initialize GBEST solution
  pso_result_t solution;
  // allocate memory for the best position buffer
  solution.gbest = malloc(settings.dim * sizeof(double));

  // run optimization algorithm
  pso_solve(obj_fun, NULL, &solution, &settings);

  // free the gbest buffer
  free(solution.gbest);

  return 0;

}
