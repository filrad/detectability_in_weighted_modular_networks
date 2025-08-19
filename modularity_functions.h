struct weighted_network{
  int N;
  int **bond;
  double **weights;
  int *in_degree;
  int *out_degree;
  int E;
  double M;
  //for poisson degree distributions
  double k_in;
  double k_out;
  //
  double w_in;
  double w_out;
  //
  //set model = 0 for unweighted, set model =1 for Poisson weights, and model = 2 for exponential weights
  int model;
};


struct spectral_modularity {
  int N;
  double *evector;
  double eval;
  double *tmp_evector;
  double err;
  double precision;
  int iterations;
  int max_iterations;
  double order_parameter;
};


void evaluate_cluster_reconstruction (struct spectral_modularity *S);
void allocate_memory_spectral (struct spectral_modularity *S);
void deallocate_memory_spectral (struct spectral_modularity *S);
void single_iteration_spectral_modularity (struct spectral_modularity *S, struct weighted_network *W);void spectral_modularity_maximization (struct spectral_modularity *S, struct weighted_network *W);
void spectral_modularity_maximization (struct spectral_modularity *S, struct weighted_network *W);
void print_spectral(struct spectral_modularity *S);




void allocate_memory (struct weighted_network *W);
void deallocate_memory (struct weighted_network *W);

int poisson_random_variable (double lambda);
double exponential_random_variable (double lambda);
int geometric_random_variable (double lambda);
int signed_random_variable (double lambda);




void generate_poisson_degree_sequences (struct weighted_network *W);
void print_network (struct weighted_network *W);

int check_connection (int **bond, int n, int m);
void generate_internal_connections(struct weighted_network *W, int i_min, int i_max);
void generate_external_connections(struct weighted_network *W, int i_min, int i_max, int j_min, int j_max);
void generate_unweighted_graph (struct weighted_network *W);
void assign_weights_to_connections (struct weighted_network *W);

void generate_binomial_network (struct weighted_network *W);
