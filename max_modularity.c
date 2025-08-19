#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "mt64.h"
#include "modularity_functions.h"





int main (int argc, char **argv)
{



  int pid_id= time(NULL) * getpid();
  init_genrand64((unsigned)pid_id);

  
  clock_t start, end;
  double cpu_time_used;

  struct weighted_network *W = (struct weighted_network*)malloc(1*sizeof(struct weighted_network));
  struct spectral_modularity *S = (struct spectral_modularity*)malloc(1*sizeof(struct spectral_modularity));

  int N = atoi(argv[1]);

  W->N = 2*N;
  S->N = W->N;
  
  S->precision = 1e-3;
  S->max_iterations = 100000;

  W->model = atoi(argv[2]);
  
  W->k_in = atof(argv[3]);
  W->k_out = atof(argv[4]);
  W->w_in = atof(argv[5]);
  W->w_out = atof(argv[6]);

 


  

  
  allocate_memory (W);
  allocate_memory_spectral (S);
  

  start = clock();
  ////////////////////////////
  //generate_poisson_degree_sequences (W);
  //generate_unweighted_graph (W);

  generate_binomial_network (W);
  
  assign_weights_to_connections (W);

  ////////////////////////////
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //printf("#generate graph model : %g\n",cpu_time_used); fflush(stdout);


  //printf("#E = %d %g\n",W->E, (double)W->E/(double)W->N);
  
  //print_network (W);


  
  start = clock();
  ////////////////////////////
  spectral_modularity_maximization (S, W);
  ////////////////////////////
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //printf("#modularity maximization : %g\n",cpu_time_used); fflush(stdout);


  evaluate_cluster_reconstruction (S);

  

  printf("%d %d %g %g %g %g\t",W->N, W->model,W->k_in,W->k_out,W->w_in,W->w_out);
  printf("%g %g %g %d\n",S->eval,S->order_parameter,S->err,S->iterations);
  
  
  //print_spectral(S);
  
  
  deallocate_memory (W);
  deallocate_memory_spectral (S);


  return 0;
}
