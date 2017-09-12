//Daniel Wilson jw563@exeter.ac.uk 27/02/17
//This program uses the Ising model to model a ferromagnetic material over a temperature range. The model uses a lattice of size D*D with periodic boundary conditions to simulate an infinte lattice
//main.c takes a command line interface directing the method of checking the lattice (ordered or random) and the lattice size
//The program stores NSAVE values of M and E and after NSAVE iterations of the metropilis algorithm the gradient of the most recent NSAVE values are checked. If below GRAD_TOL, the set is accepted and the temperature is incramented forwards.
//Upon reaching T_MAX, the program stores the value of M, Chi and Cap at T_C for use in finding the critical exponents of the system.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>

#define RANDOM_MAX 0x7FFFFFFF
#define S(x,y) S[(x+D)%D][(y+D)%D]
#define T_MAX 5
#define T_MIN 0
#define mu0B 0
#define T_C 2.5
#define k_B 1.38e-23
#define ITERATIONS 100
#define EVOLVE_MAX 1000000
#define NSAVE 10000
#define GRAD_TOL 1e-7
#define J k_B*T_C*0.4406867935 // J/K
#define T_STEP ((double)T_MAX-(double)T_MIN)/ITERATIONS

#define GNUPLOT_DATA "output_data.txt"

typedef struct properties {
double M_save[NSAVE];
double E_save[NSAVE];
double S_M;
double S_E;
double chi_max;
double cap_max;
double chi_T_C;
double cap_T_C;
double M_T_C;
} Properties;

int D;

//allocates memory for the struct properties, using calloc() to allocate cleared memory, and exiting the program if there is insufficient memory.
static void *xmalloc(size_t n) {
  void *p = calloc(1,n);
  if (p == NULL) {
    fprintf(stderr, "Out of memory!\n");
    exit(1);
  }
  return p;
}

//evaluates the properties of a set of NSAVE values for M and E which are considered to be randomly fluctuating about an average. The function sums the squares of E and M, and uses the sums from grad_check(), for use in calculating chi (susceptibility) and cap (heat capacity). The variable k is Onsager's exact solution, and is plotted alongside M. The maximum values for chi and cap, as well as the value for M at T_C are also stored.
static void prop_eval(Properties *props, double T, FILE *foutput) {

  double S_MM=0.0,S_EE=0.0;
  double chi, cap, k;
  for (int i=0; i<NSAVE; i++) {
      S_MM += (props->M_save[i]*props->M_save[i]) / NSAVE;
      S_EE += (props->E_save[i]*props->E_save[i]) / NSAVE;
  }

  props->S_M = props->S_M/ NSAVE;
  props->S_E = props->S_E / NSAVE;

  chi = (1/T) * (S_MM - (props->S_M*props->S_M)) * (D*D);
  if (chi > props->chi_max) {
    props->chi_max = chi;
    props->chi_T_C = T;
  }

  cap = (1/(T*T*J*J)) * (S_EE - (props->S_E*props->S_E)) * (D*D);
  if (cap > props->cap_max) {
    props->cap_max = cap;
    props->cap_T_C = T;
  }

  k = 1/(sinh((2*J)/(k_B*T)));
  if (k > 1) k = 0;
  else k = pow((1-k*k*k*k),0.125);

  fprintf(foutput, "%g\t%g\t%g\t%g\t%g\n", T , props->S_M, k,chi,cap );

  if (fabs(T-T_C) < T_STEP/2) props->M_T_C = props->S_M;

}

//checks the gradient of the values within the arrays M_save and E_save. If the gradients are both below the threshold GRAD_TOL, then the set of data is accepted
static int grad_check(Properties *props, double T,int i) {

  if (i>=(NSAVE-1)) {

    double M_grad, E_grad;
    double S_xx = 0, S_Mx = 0, S_x = 0, S_Ex = 0;
    int x;

    props->S_M = 0.0;
    props->S_E = 0.0;

    for (int j =0; j<NSAVE; j++) {
      x = (1+i+j)%NSAVE;
      S_x += x;
      S_xx += x*x;
      S_Mx += x*fabs(props->M_save[x]);
      S_Ex += x*props->E_save[x];

      props->S_M += fabs(props->M_save[x]);
      props->S_E += props->E_save[x];
    }

    M_grad = fabs((NSAVE*S_Mx - S_x*props->S_M) / (NSAVE*S_xx - S_x*S_x));
    E_grad = fabs((NSAVE*S_Ex - S_x*props->S_E) / (NSAVE*S_xx - S_x*S_x));

    if ((M_grad < GRAD_TOL) && (E_grad < GRAD_TOL)) {
      printf("T = %5gK\ti = %7d\r",T,i+1);
      fflush(stdout);
      return 1;
    }
  }
return 0;
}

//calculates and stores the values of E and M from the pervious evolution in the oldest element of the array M_save/E_save
static void prop_store(Properties *props, short **S, int i,double T) {
  double M_total = 0.0;
  double E_total = 0.0;

  for (int x=0; x<D; x++) {
    for (int y=0; y<D; y++) {
      M_total += S(x,y);
      E_total -= S(x,y)*( J*(S(x+1,y)+S(x,y+1)) + mu0B);
    }
  }
  props->M_save[i%NSAVE] = M_total / (D*D);
  props->E_save[i%NSAVE] = E_total / (D*D);
}

//invokes the metropilis algorithm for deciding when/if to filp the spin of a specific site
static void metro(Properties *props, short **S, double T, int x, int y) {
  double Delta_E;
  double U;
  Delta_E = S(x,y)*( 2.0*J*(S(x-1,y)+S(x,y-1)+S(x+1,y)+S(x,y+1)) + mu0B );

  if (Delta_E <= 0.0) {
    S(x,y) = -S(x,y);
  }
  else {
    U = random()/(double)RANDOM_MAX;
    if (U < exp(-(double)Delta_E / ((double)k_B * T)) ) {
      S(x,y) = -S(x,y);
    }
  }
}

//evolves the grid either randomly ('rnd') or in order ('ord'). After D*D metropilis evolutions, the properties of the grid are stored, and then the gradient of the variables is calculated to see if they are accepted. If so, the function breaks and the next temperature step is added in main(). Otherwise the grid continues to evolve until EVOLVE_MAX is reached
static void evolve_grid(short **S, Properties *props, double T,char option) {
  int UX, UY;
  switch(option) {

    case 'o' :
    for(int i=0; i<EVOLVE_MAX; i++) {
      for (int x=0; x<D; x++) {
        for (int y=0; y<D; y++) {
          metro(props,S,T,x,y);
        }
      }
      prop_store(props,S,i,T);
      if (grad_check(props,T,i) == 1) break;
    }
    break;

    case 'r' :
    for(int i=0; i<EVOLVE_MAX; i++) {
      for (int j=0; j< (D*D); j++){
        UX = random()%D;
        UY = random()%D;
        metro(props,S,T,UX,UY);
      }
      prop_store(props,S,i,T);
      if (grad_check(props,T,i) == 1) break;
    }
    break;

    default :
    printf("Unknown input.\nPlease input one of the following:\n'r' (Random selection of sites)\n'o' (Ordered selection of sites)\nFollowed by the size of the grid.\ne.g. ./program o 20\n");
    exit(1);
  }
}

//randomly allocates each lattice point a spin of +1.0 (UP) or -1.0 (DOWN)
static void set_grid(short **S, Properties *props) {
  int U;
  for (int x=0; x<D; x++) {
    for(int y=0; y<D; y++) {
      U = random()%2;
      if (U==1) S(x,y) = 1.0;
      else S(x,y) = -1.0;
    }
  }
}

//checks for the correct number of inputs from the command line and for a valid grid size. The while loop ensures that the first evolution results in a spin of +1.0 - due to the random nature of the grid allocation, and the lack of an external field, the system can sometimes find equilibrium at M != 1, which we are not interested in. The temperature is then stepped forwards up to T_MAX, the relevant data is added to 'criticaldata.txt' for calculating the critical exponents, and the graphs for M, chi and cap are plotted.
int main (int argc, const char *argv[]) {

  if (argc!=3) {
    printf("Unknown input.\nPlease input one of the following:\n'r' (Random selection of sites)\n'o' (Ordered selection of sites)\nFollowed by the size of the grid.\ne.g. ./program o 20\n");
    exit(1);
  }

  char option = argv[1][0];

  D = atoi (argv[2]);

  Properties *props = xmalloc(sizeof(Properties));
  srandom((unsigned int)time(NULL));

  short **S;
  S = xmalloc(D*sizeof(*S));
  for (int i=0;i<D;i++) {
    S[i] = xmalloc(D*sizeof(*S[i]));
  }
  //short S[D][D];

  while (1) {
    set_grid(S,props);
    evolve_grid(S,props, T_MIN + T_STEP ,option);
    if (!(props->S_M/NSAVE < 0.9)) break;
  }

  double T = T_MIN;

  FILE *foutput = fopen (GNUPLOT_DATA, "w");

  for (int i = 0; i < ITERATIONS; i++) {
    T += T_STEP;
    evolve_grid(S,props,T,option);
    prop_eval(props,T,foutput);
  }

  fclose(foutput);

  FILE *criticaldata = fopen("criticaldata.txt", "a");
  fprintf(criticaldata, "%d\t%g\t%g\t%g\t%g\t%g\t%g\n", D, log((double)D), log(props->M_T_C), log(props->chi_max), log(props->cap_max), props->chi_T_C, props->cap_T_C);
  fclose(criticaldata);

  free(props);

  for (int i=0;i<D;i++) {
    free(S[i]);
  }
  free(S);


  return 0;

}
