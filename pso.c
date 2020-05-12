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


#include <stdlib.h> // for rand() stuff
#include <stdio.h> // for printf
#include <time.h> // for time()
#include <math.h> // for cos(), pow(), sqrt() etc.
#include <float.h> // for DBL_MAX
#include <string.h> // for mem*

#include "pso.h"

// generates a double between (0, 1)
#define RNG_UNIFORM() (rand()/(double)RAND_MAX)

// generate an int between 0 and s (exclusive)
#define RNG_UNIFORM_INT(s) (rand()%s)

// function type for the different inform functions
typedef void (*inform_fun_t)(int *comm, double **pos_nb,
                             double **pos_b, double *fit_b,
                             double *gbest, int improved,
                             pso_settings_t *settings);

// function type for the different inertia calculation functions
typedef double (*inertia_fun_t)(int step, pso_settings_t *settings);


//==============================================================
// calulate swarm size based on dimensionality
int pso_calc_swarm_size(int dim) {
    int size = 10. + 2. * sqrt(dim);
    return (size > PSO_MAX_SIZE ? PSO_MAX_SIZE : size);
}


//==============================================================
//          INERTIA WEIGHT UPDATE STRATEGIES
//==============================================================
// calculate linearly decreasing inertia weight
double calc_inertia_lin_dec(int step, pso_settings_t *settings) {

    int dec_stage = 3 * settings->steps / 4;
    if (step <= dec_stage)
        return settings->w_min + (settings->w_max - settings->w_min) *	\
            (dec_stage - step) / dec_stage;
    else
        return settings->w_min;
}



//==============================================================
//          NEIGHBORHOOD (COMM) MATRIX STRATEGIES
//==============================================================
// global neighborhood
void inform_global(int *comm, double **pos_nb,
		   double **pos_b, double *fit_b,
		   double *gbest, int improved,
		   pso_settings_t *settings)
{

    int i;
    // all particles have the same attractor (gbest)
    // copy the contents of gbest to pos_nb
    for (i=0; i<settings->size; i++)
        memmove((void *)pos_nb[i], (void *)gbest,
                sizeof(double) * settings->dim);

}


// ===============================================================
// general inform function :: according to the connectivity
// matrix COMM, it copies the best position (from pos_b) of the
// informers of each particle to the pos_nb matrix
void inform(int *comm, double **pos_nb, double **pos_b, double *fit_b,
	    int improved, pso_settings_t * settings)
{
    int i, j;
    int b_n; // best neighbor in terms of fitness

    // for each particle
    for (j=0; j<settings->size; j++) {
        b_n = j; // self is best
        // who is the best informer??
        for (i=0; i<settings->size; i++)
            // the i^th particle informs the j^th particle
            if (comm[i*settings->size + j] && fit_b[i] < fit_b[b_n])
                // found a better informer for j^th particle
                b_n = i;
        // copy pos_b of b_n^th particle to pos_nb[j]
        memmove((void *)pos_nb[j],
                (void *)pos_b[b_n],
                sizeof(double) * settings->dim);
    }
}




// =============
// ring topology
// =============

// topology initialization :: this is a static (i.e. fixed) topology
void init_comm_ring(int *comm, pso_settings_t * settings) {
    int i;
    // reset array
    memset((void *)comm, 0, sizeof(int)*settings->size*settings->size);

    // choose informers
    for (i=0; i<settings->size; i++) {
        // set diagonal to 1
        comm[i*settings->size+i] = 1;
        if (i==0) {
            // look right
            comm[i*settings->size+i+1] = 1;
            // look left
            comm[(i+1)*settings->size-1] = 1;
        } else if (i==settings->size-1) {
            // look right
            comm[i*settings->size] = 1;
            // look left
            comm[i*settings->size+i-1] = 1;
        } else {
            // look right
            comm[i*settings->size+i+1] = 1;
            // look left
            comm[i*settings->size+i-1] = 1;
        }

    }

}




void inform_ring(int *comm, double **pos_nb,
		 double **pos_b, double *fit_b,
		 double *gbest, int improved,
		 pso_settings_t * settings)
{

    // update pos_nb matrix
    inform(comm, pos_nb, pos_b, fit_b, improved, settings);

}

// ============================
// random neighborhood topology
// ============================
void init_comm_random(int *comm, pso_settings_t * settings) {

    int i, j, k;
    // reset array
    memset((void *)comm, 0, sizeof(int)*settings->size*settings->size);

    // choose informers
    for (i=0; i<settings->size; i++) {
        // each particle informs itself
        comm[i*settings->size + i] = 1;
        // choose kappa (on average) informers for each particle
        for (k=0; k<settings->nhood_size; k++) {
            // generate a random index
            j = RNG_UNIFORM_INT(settings->size);
            // particle i informs particle j
            comm[i*settings->size + j] = 1;
        }
    }
}



void inform_random(int *comm, double **pos_nb,
		   double **pos_b, double *fit_b,
		   double *gbest, int improved,
		   pso_settings_t * settings)
{


    // regenerate connectivity??
    if (!improved)
        init_comm_random(comm, settings);
    inform(comm, pos_nb, pos_b, fit_b, improved, settings);

}




//==============================================================
// create pso settings
pso_settings_t *pso_settings_new(int dim, double range_lo, double range_hi) {
    pso_settings_t *settings = (pso_settings_t *)malloc(sizeof(pso_settings_t));
    if (settings == NULL) { return NULL; }

    // set some default values
    settings->dim = dim;
    settings->goal = 1e-5;

    // set up the range arrays
    settings->range_lo = (double *)malloc(settings->dim * sizeof(double));
    if (settings->range_lo == NULL) { free(settings); return NULL; }

    settings->range_hi = (double *)malloc(settings->dim * sizeof(double));
    if (settings->range_hi == NULL) { free(settings); free(settings->range_lo); return NULL; }

    for (int i=0; i<settings->dim; i++) {
        settings->range_lo[i] = range_lo;
        settings->range_hi[i] = range_hi;
    }

    settings->size = pso_calc_swarm_size(settings->dim);
    settings->print_every = 1000;
    settings->steps = 100000;
    settings->c1 = 1.496;
    settings->c2 = 1.496;
    settings->w_max = PSO_INERTIA;
    settings->w_min = 0.3;

    settings->clamp_pos = 1;
    settings->nhood_strategy = PSO_NHOOD_RING;
    settings->nhood_size = 5;
    settings->w_strategy = PSO_W_LIN_DEC;

    return settings;
}

// destroy PSO settings
void pso_settings_free(pso_settings_t *settings) {
    free(settings->range_lo);
    free(settings->range_hi);
    free(settings);
}


double **pso_matrix_new(int size, int dim) {
    double **m = (double **)malloc(size * sizeof(double *));
    for (int i=0; i<size; i++) {
        m[i] = (double *)malloc(dim * sizeof(double));
    }
    return m;
}

void pso_matrix_free(double **m, int size) {
    for (int i=0; i<size; i++) {
        free(m[i]);
    }
    free(m);
}


//==============================================================
//                     PSO ALGORITHM
//==============================================================
void pso_solve(pso_obj_fun_t obj_fun, void *obj_fun_params,
	       pso_result_t *solution, pso_settings_t *settings)
{
    // Particles
    double **pos = pso_matrix_new(settings->size, settings->dim); // position matrix
    double **vel = pso_matrix_new(settings->size, settings->dim); // velocity matrix
    double **pos_b = pso_matrix_new(settings->size, settings->dim); // best position matrix
    double *fit = (double *)malloc(settings->size * sizeof(double));
    double *fit_b = (double *)malloc(settings->size * sizeof(double));
    // Swarm
    double **pos_nb = pso_matrix_new(settings->size, settings->dim); // what is best informed
    // position for each particle
    int *comm = (int *)malloc(settings->size * settings->size * sizeof(int));
    // rows : those who inform
    // cols : those who are informed
    int improved = 0; // whether solution->error was improved during
    // the last iteration

    int i, d, step;
    double a, b; // for matrix initialization
    double rho1, rho2; // random numbers (coefficients)
    // initialize omega using standard value
    double w = PSO_INERTIA;
    inform_fun_t inform_fun = NULL; // neighborhood update function
    inertia_fun_t calc_inertia_fun = NULL; // inertia weight update function

    // initialize random seed
    srand(time(NULL));

    // SELECT APPROPRIATE NHOOD UPDATE FUNCTION
    switch (settings->nhood_strategy)
        {
        case PSO_NHOOD_GLOBAL:
            // comm matrix not used
            inform_fun = inform_global;
            break;
        case PSO_NHOOD_RING:
            init_comm_ring(comm, settings);
            inform_fun = inform_ring;
            break;
        case PSO_NHOOD_RANDOM:
            init_comm_random(comm, settings);
            inform_fun = inform_random;
            break;
        default:
            // use global as the default
            inform_fun = inform_global;
            break;
        }

    // SELECT APPROPRIATE INERTIA WEIGHT UPDATE FUNCTION
    switch (settings->w_strategy)
        {
            /* case PSO_W_CONST : */
            /*     calc_inertia_fun = calc_inertia_const; */
            /*     break; */
        case PSO_W_LIN_DEC :
            calc_inertia_fun = calc_inertia_lin_dec;
            break;
        }

    // INITIALIZE SOLUTION
    solution->error = DBL_MAX;

    // SWARM INITIALIZATION
    // for each particle
    for (i=0; i<settings->size; i++) {
        // for each dimension
        for (d=0; d<settings->dim; d++) {
            // generate two numbers within the specified range
            a = settings->range_lo[d] + (settings->range_hi[d] - settings->range_lo[d]) * \
                RNG_UNIFORM();
            b = settings->range_lo[d] + (settings->range_hi[d] - settings->range_lo[d]) *	\
                RNG_UNIFORM();
            // initialize position
            pos[i][d] = a;
            // best position is the same
            pos_b[i][d] = a;
            // initialize velocity
            vel[i][d] = (a-b) / 2.;
        }
        // update particle fitness
        fit[i] = obj_fun(pos[i], settings->dim, obj_fun_params);
        fit_b[i] = fit[i]; // this is also the personal best
        // update gbest??
        if (fit[i] < solution->error) {
            // update best fitness
            solution->error = fit[i];
            // copy particle pos to gbest vector
            memmove((void *)solution->gbest, (void *)pos[i],
                    sizeof(double) * settings->dim);
        }

    }

    // RUN ALGORITHM
    for (step=0; step<settings->steps; step++) {
        // update current step
        settings->step = step;
        // update inertia weight
        // do not bother with calling a calc_w_const function
        if (calc_inertia_fun != NULL) {
            w = calc_inertia_fun(step, settings);
        }
        // check optimization goal
        if (solution->error <= settings->goal) {
            // SOLVED!!
            if (settings->print_every)
                printf("Goal achieved @ step %d (error=%.3e) :-)\n", step, solution->error);
            break;
        }

        // update pos_nb matrix (find best of neighborhood for all particles)
        inform_fun(comm, (double **)pos_nb, (double **)pos_b,
                   fit_b, solution->gbest, improved, settings);
        // the value of improved was just used; reset it
        improved = 0;

        // update all particles
        for (i=0; i<settings->size; i++) {
            // for each dimension
            for (d=0; d<settings->dim; d++) {
                // calculate stochastic coefficients
                rho1 = settings->c1 * RNG_UNIFORM();
                rho2 = settings->c2 * RNG_UNIFORM();
                // update velocity
                vel[i][d] = w * vel[i][d] +	\
                    rho1 * (pos_b[i][d] - pos[i][d]) +	\
                    rho2 * (pos_nb[i][d] - pos[i][d]);
                // update position
                pos[i][d] += vel[i][d];
                // clamp position within bounds?
                if (settings->clamp_pos) {
                    if (pos[i][d] < settings->range_lo[d]) {
                        pos[i][d] = settings->range_lo[d];
                        vel[i][d] = 0;
                    } else if (pos[i][d] > settings->range_hi[d]) {
                        pos[i][d] = settings->range_hi[d];
                        vel[i][d] = 0;
                    }
                } else {
                    // enforce periodic boundary conditions
                    if (pos[i][d] < settings->range_lo[d]) {

                        pos[i][d] = settings->range_hi[d] - fmod(settings->range_lo[d] - pos[i][d],
                                                                 settings->range_hi[d] - settings->range_lo[d]);
                        vel[i][d] = 0;

                    } else if (pos[i][d] > settings->range_hi[d]) {

                        pos[i][d] = settings->range_lo[d] + fmod(pos[i][d] - settings->range_hi[d],
                                                                 settings->range_hi[d] - settings->range_lo[d]);
                        vel[i][d] = 0;
                    }
                }

            }

            // update particle fitness
            fit[i] = obj_fun(pos[i], settings->dim, obj_fun_params);
            // update personal best position?
            if (fit[i] < fit_b[i]) {
                fit_b[i] = fit[i];
                // copy contents of pos[i] to pos_b[i]
                memmove((void *)pos_b[i], (void *)pos[i],
                        sizeof(double) * settings->dim);
            }
            // update gbest??
            if (fit[i] < solution->error) {
                improved = 1;
                // update best fitness
                solution->error = fit[i];
                // copy particle pos to gbest vector
                memmove((void *)solution->gbest, (void *)pos[i],
                        sizeof(double) * settings->dim);
            }
        }

        if (settings->print_every && (step % settings->print_every == 0))
            printf("Step %d (w=%.2f) :: min err=%.5e\n", step, w, solution->error);

    }

    // free resources
    pso_matrix_free(pos, settings->size);
    pso_matrix_free(vel, settings->size);
    pso_matrix_free(pos_b, settings->size);
    pso_matrix_free(pos_nb, settings->size);
    free(comm);
    free(fit);
    free(fit_b);
}
