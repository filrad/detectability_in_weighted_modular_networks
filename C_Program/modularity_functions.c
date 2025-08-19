#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "mt64.h"
#include "modularity_functions.h"






void evaluate_cluster_reconstruction (struct spectral_modularity *S)
{

  int n, N;
  double cluster1, cluster2;
  cluster1 = cluster2 = 0.0;

  N = S->N/2;
  
  for(n=1;n<=N;n++) cluster1 += S->evector[n];
  for(n=N+1;n<=2*N;n++) cluster2 += S->evector[n];

  S->order_parameter = fabs(cluster1) + fabs(cluster2);
  S->order_parameter = S->order_parameter / sqrt((double)S->N);

  
}



///////////////////////////////////
void print_spectral(struct spectral_modularity *S)
{
  int n;
  printf("# %g\n",S->eval);
  for(n=1;n<=S->N;n++) printf("%d %g\n",n,S->evector[n]);
}
///////////////////////////////////


////////////////////////////////////
void spectral_modularity_maximization (struct spectral_modularity *S, struct weighted_network *W)
{

  S->iterations = 0.0;
  
  int n;
  for (n=1;n<=S->N;n++) S->evector[n] = genrand64_real3();
  S->eval = 0.0;
  for(n=1;n<=S->N;n++) S->eval += S->evector[n]*S->evector[n];
  S->eval = sqrt(S->eval);
  for(n=1;n<=S->N;n++) S->evector[n] = S->evector[n] / S->eval;

  S->err = 1.0;
  while ((S->err > S->precision) && (S->iterations < S->max_iterations) )
    {
      single_iteration_spectral_modularity (S, W);
      //printf("#%d %g\n",S->iterations, S->err);
    }
 
}













//////////////////////////////////////
void single_iteration_spectral_modularity (struct spectral_modularity *S, struct weighted_network *W)
{

  int n, m, i;
  double err, summ, sign;

  S->iterations += 1;

  summ = 0.0;
  for(n=1;n<=S->N;n++)  if(W->bond[0][n]>0) summ += (S->evector[n] * W->weights[0][n]);
  summ /= W->M;
  
  //printf("#W = %g\n",W->M);
  //printf("#S = %g\n",summ);
  
  for(n=1;n<=S->N;n++)
    {
      if(W->bond[0][n] > 0){
	S->tmp_evector[n] = - summ * W->weights[0][n];
	for(i=1;i<=W->bond[0][n];i++)
	  {
	    m = W->bond[n][i];
	    S->tmp_evector[n] += S->evector[m] * W->weights[n][i];
	  }
      }
    }


  //flip sign for convergence
  i = 1;
  while(W->bond[0][i]==0) i++;
  sign = 1.0;
  if (S->tmp_evector[i] < 0.0) sign= -1.0;
  //printf("#flip %d %g\n",i, sign);
  for(n=1;n<=S->N;n++) S->tmp_evector[n] *= sign;
  
  
  S->eval = 0.0;
  for(n=1;n<=S->N;n++) if(W->bond[0][n]>0) S->eval += S->tmp_evector[n]*S->tmp_evector[n];
  S->eval = sqrt(S->eval);




  
  S->err = 0.0;
   for(n=1;n<=S->N;n++)
     {
       if(W->bond[0][n]>0){
	 S->tmp_evector[n] = S->tmp_evector[n] / S->eval;
	 err = fabs(S->tmp_evector[n]-S->evector[n]);
	 if (err > S->err) S->err = err;
	 S->evector[n] = S->tmp_evector[n];
       }
     }
}



/////////////////
void allocate_memory_spectral (struct spectral_modularity *S)
{  
  S->evector = (double *)malloc((S->N+1)*sizeof(double *));
  S->tmp_evector = (double *)malloc((S->N+1)*sizeof(double *));
}

void deallocate_memory_spectral (struct spectral_modularity *S)
{  
  free(S->evector);
  free(S->tmp_evector);
  free(S);
}
//////////////////

/////////////////
void allocate_memory (struct weighted_network *W)
{

  int i;
  W->bond = (int **)malloc((W->N+1)*sizeof(int **));
  W->weights = (double **)malloc((W->N+1)*sizeof(double **));
  W->bond[0] = (int *)malloc((W->N+1)*sizeof(int *));
  W->weights[0] = (double *)malloc((W->N+1)*sizeof(double *));
  for(i=1;i<=W->N;i++)
    {
      W->bond[0][i] = 0;
      W->weights[0][i] = 0.0;
      W->bond[i] = (int *)malloc(1*sizeof(int));
      W->weights[i] = (double *)malloc(1*sizeof(double));
    }
  W->in_degree = (int *)malloc((W->N+1)*sizeof(int *));
  W->out_degree = (int *)malloc((W->N+1)*sizeof(int *));
}
//////////////////
void deallocate_memory (struct weighted_network *W)
{
  int i;
  for(i=0;i<=W->N;i++)
    {
      free(W->bond[i]);
      free(W->weights[i]);
    }

  free(W->bond);
  free(W->weights);
  free(W->in_degree);
  free(W->out_degree);
  free(W);
}
//////////////////

/*
int poisson_random_variable (double lambda) {
  int n = 0;
  double limit; 
  double x; 
  limit = exp(-lambda);
  x = genrand64_real3(); 
  while (x > limit) {
    n++;
    x *= genrand64_real3(); 
  }
  return n;
}
*/


	
//Poisson function -- returns a single Poisson random variable
int poisson_random_variable (double lambda)
{
  double exp_lambda = exp(-lambda); //constant for terminating loop
  double randUni; //uniform variable
  double prodUni; //product of uniform variables
  int randPoisson; //Poisson variable
 
  //initialize variables
  randPoisson = -1;
  prodUni = 1;
  do
    {
      randUni = genrand64_real3(); //generate uniform variable
      prodUni = prodUni * randUni; //update product
      randPoisson++; //increase Poisson variable
 
    } while (prodUni > exp_lambda);
  return randPoisson;
}



double exponential_random_variable (double lambda){
    double u;
    u = genrand64_real3();
    return -log(1- u) / lambda;
}
//////////////////



int geometric_random_variable (double lambda){
  double u, tmp;
  int a;
  if (lambda <=0 ) return -1;
  if (lambda == 1 ) return 1;

  u = log(genrand64_real3());
  tmp = log(1.0-lambda);
  a = 1 + u/tmp;
  return a;
 }

//////////////////


int signed_random_variable (double lambda){
    double u;
    u = genrand64_real3();
    if (u < lambda) return 1;
    return -1;
}





void print_network (struct weighted_network *W)
{
  int i, j;
  //for(i=1;i<=W->N;i++) printf("%d %d %d\n", i, W->in_degree[i], W->out_degree[i]);
  for(i=1;i<=W->N;i++)
    {
      printf("%d",i);
      for(j=1;j<=W->bond[0][i];j++) printf(" %d(%g) ",W->bond[i][j],W->weights[i][j]);
      printf("\n");
    }
  
}
//////////////////
void generate_poisson_degree_sequences (struct weighted_network *W)
{
  int i;
  for(i=1;i<=W->N;i++)
    {
      W->in_degree[i] =  poisson_random_variable (W->k_in);
      W->out_degree[i] =  poisson_random_variable (W->k_out);
    }
}





int check_connection (int **bond, int n, int m)
{
  int i;
  for(i=1;i<=bond[0][n];i++) if (bond[n][i] == m) return i;
  return -1;
}


void generate_internal_connections(struct weighted_network *W, int i_min, int i_max)
{

  int i, j, tmp, n, m, E = 0;
  for (i=i_min;i<=i_max;i++)
    {
      E += W->in_degree[i];
    }

  int *edges = (int *)malloc((E+1)*sizeof(int));
  for(i=1;i<=E;i++) edges[i] = -1;

  
  
  E = 0;
  for (i=i_min;i<=i_max;i++)
    {
      for(j=1;j<=W->in_degree[i];j++)
	{
	  E++;
	  edges[E] = i;
	}
    }
  //random shuffle
  /*
  i = 0;
  while(i<10*E)
    {
      n = (int)((double)E*genrand64_real3())+1;
      if (n>E) n=1;
      m = (int)((double)E*genrand64_real3())+1;
      if (m>E) m=1;
      tmp = edges[n];
      edges[n] = edges[m];
      edges[m] = tmp;
      i++;
    }
  */

  for(n=1;n<=E;n++)
    {
       m = (int)((double)E*genrand64_real3())+1;
       if (m>E) m=1;
       tmp = edges[n];
       edges[n] = edges[m];
       edges[m] = tmp;
    }


  
  

 
  ////////////////

  int count = 0;
  //create connections
  i = 0;
   while(i<E-1)
     {
       i++;
       n = edges[i];
       i++;
       m = edges[i];

       if(n!=m)
	 {
	   if(check_connection (W->bond, n, m) <0){
	     //printf("%d %d %d\n",i,n,m);
	     W->bond[0][n] += 1;
	     W->bond[n] = (int *)realloc(W->bond[n], (W->bond[0][n]+1)*sizeof(int));
	     W->bond[n][W->bond[0][n]] = m;
	     W->bond[0][m] += 1;
	     W->bond[m] = (int *)realloc(W->bond[m], (W->bond[0][m]+1)*sizeof(int));
	     W->bond[m][W->bond[0][m]] = n;
	     W->E += 2;
	   }
	 }
       else{
	 count +=1;
	 //printf("%d %d\n",n,m);
       }
     }
   //printf ("#int C = %d\n",count);
  ////////////////////

  
  

  free(edges);
}









void generate_external_connections(struct weighted_network *W, int i_min, int i_max, int j_min, int j_max)
{

  int i, j, tmp, n, m, E1, E2, E;

  E1 = 0;
  for (i=i_min;i<=i_max;i++)
    {
      E1 += W->out_degree[i];
    }
  E2 = 0;
  for (j=j_min;j<=j_max;j++)
    {
      E2 += W->out_degree[j];
    }
  
  
  int *edges1 = (int *)malloc((E1+1)*sizeof(int));
  int *edges2 = (int *)malloc((E2+1)*sizeof(int));

  E1 = 0;
  for (i=i_min;i<=i_max;i++)
    {
      for(j=1;j<=W->out_degree[i];j++)
	{
	  E1++;
	  edges1[E1] = i;
	}
    }
  E2= 0;
  for (j=j_min;j<=j_max;j++)
    {
      for(i=1;i<=W->out_degree[j];i++)
	{
	  E2++;
	  edges2[E2] = j;
	}
    }


  //random shuffle
  /*
  i = 0;
  while(i<10*E1)
    {
      n = (int)((double)E1*genrand64_real3())+1;
      if (n>E1) n=1;
      m = (int)((double)E1*genrand64_real3())+1;
      if (m>E1) m=1;
      tmp = edges1[n];
      edges1[n] = edges1[m];
      edges1[m] = tmp;
      i++;
    }
  */
  for(n=1;n<=E1;n++)
    {
       m = (int)((double)E1*genrand64_real3())+1;
       if (m>E1) m=1;
       tmp = edges1[n];
       edges1[n] = edges1[m];
       edges1[m] = tmp;
    }
  
  //random shuffle
  /*
  i = 0;
  while(i<10*E2)
    {
      n = (int)((double)E2*genrand64_real3())+1;
      if (n>E2) n=1;
      m = (int)((double)E2*genrand64_real3())+1;
      if (m>E2) m=1;
      tmp = edges2[n];
      edges2[n] = edges2[m];
      edges2[m] = tmp;
      i++;
    }
  */
  for(n=1;n<=E2;n++)
    {
       m = (int)((double)E2*genrand64_real3())+1;
       if (m>E2) m=1;
       tmp = edges2[n];
       edges2[n] = edges2[m];
       edges2[m] = tmp;
    }
  ////////////////


  int count=0;

  //create connections
  i = 0;
  E = E1;
  if(E2<E) E = E2;
   while(i<E-1)
     {
       i++;
       n = edges1[i];
       m = edges2[i];

       if(n!=m)
	 {
	   if(check_connection (W->bond, n, m) <0){
	     //printf("%d %d %d\n",i,n,m);
	     W->bond[0][n] += 1;
	     W->bond[n] = (int *)realloc(W->bond[n], (W->bond[0][n]+1)*sizeof(int));
	     W->bond[n][W->bond[0][n]] = m;
	     
	     W->bond[0][m] += 1;
	     W->bond[m] = (int *)realloc(W->bond[m], (W->bond[0][m]+1)*sizeof(int));
	     W->bond[m][W->bond[0][m]] = n;
	     
	     W->E += 2;
	   }
	 }
       else{
	 count+=1;
       }
     }
  ////////////////////

   //printf("# out C = %d\n",count);
  
  

  free(edges1);
  free(edges2);
}


////////////////////////////////
void generate_unweighted_graph (struct weighted_network *W)
{
  int N = W->N/2;
  int var = 0;
  W->E = 0;
  generate_internal_connections(W, 1, N);
  //printf("# group 1 = %d\n",W->E);
  var = W->E;
  generate_internal_connections(W, N+1, 2*N);
  //printf("# group 2 = %d\n",W->E-var);
  var = W->E;
  generate_external_connections(W, 1, N, N+1, 2*N);
  //printf("# groups 1,2 = %d\n",W->E-var);
}


///////////////////////////


void assign_weights_to_connections (struct weighted_network *W)
{
  int N, n, m, i, j;
  double w, lambda;
  N = W->N/2;
  W->M = 0.0;
  for(n=1;n<=W->N;n++) W->weights[n] = (double *)realloc(W->weights[n], (W->bond[0][n]+1)*sizeof(double));
  for(n=1;n<=W->N;n++)
    {
      for(i=1;i<=W->bond[0][n];i++)
	{
	  m = W->bond[n][i];
	  if (n<m)
	    {
	      lambda = W->w_in;
	      if ((n<=N && m>N) || (n>N && m<=N)) lambda = W->w_out;
	      w = lambda; //dirac
	      if(W->model==1) w = (double) poisson_random_variable (lambda); //poisson
	      if(W->model==2) w = exponential_random_variable (1.0/lambda); //exponential
	      if(W->model==3) w = (double) geometric_random_variable (1.0/lambda); //geometric
	      if(W->model==4) w = (double) signed_random_variable (0.5*(lambda+1.0)); //sign   
	      
	      W->weights[n][i] = w;
	      W->weights[m][check_connection (W->bond, m, n)] = w;
	      
	      W->weights[0][n] += w;
	      W->weights[0][m] += w;
	      W->M += 2.0*w;
	    }
	}
    }
}






void generate_binomial_network (struct weighted_network *W)
{

  int n, m, N = W->N/2;
  double pin = W->k_in / (double)(N-1);
  double pout = W->k_out / (double)(N);

  W->E = 0;
  //internal 1
  for (n=1;n<N;n++)
    {
      for(m=n+1;m<=N;m++)
	{
	  if(genrand64_real3()<pin)
	    {
	      W->bond[0][n] += 1;
	      W->bond[n] = (int *)realloc(W->bond[n], (W->bond[0][n]+1)*sizeof(int));
	      W->bond[n][W->bond[0][n]] = m;
	      
	      W->bond[0][m] += 1;
	      W->bond[m] = (int *)realloc(W->bond[m], (W->bond[0][m]+1)*sizeof(int));
	      W->bond[m][W->bond[0][m]] = n;
	     
	      W->E += 2;
	    }
	}
    }

  //internal 2
  for (n=N+1;n<2*N;n++)
    {
      for(m=n+1;m<=2*N;m++)
	{
	  if(genrand64_real3()<pin)
	    {
	      W->bond[0][n] += 1;
	      W->bond[n] = (int *)realloc(W->bond[n], (W->bond[0][n]+1)*sizeof(int));
	      W->bond[n][W->bond[0][n]] = m;
	      
	      W->bond[0][m] += 1;
	      W->bond[m] = (int *)realloc(W->bond[m], (W->bond[0][m]+1)*sizeof(int));
	      W->bond[m][W->bond[0][m]] = n;
	     
	      W->E += 2;
	    }
	}
    }


  //external
  for (n=1;n<=N;n++)
    {
      for(m=N+1;m<=2*N;m++)
	{
	  if(genrand64_real3()<pout)
	    {
	      W->bond[0][n] += 1;
	      W->bond[n] = (int *)realloc(W->bond[n], (W->bond[0][n]+1)*sizeof(int));
	      W->bond[n][W->bond[0][n]] = m;
	      
	      W->bond[0][m] += 1;
	      W->bond[m] = (int *)realloc(W->bond[m], (W->bond[0][m]+1)*sizeof(int));
	      W->bond[m][W->bond[0][m]] = n;
	     
	      W->E += 2;
	    }
	}
    }

  
 
  
}
