#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mersenne.h"
/*use - Wall when compioling*/
/* Implement integrated autocorrelation*/
/*Measure magnetization*/
/* s=exp(itheta)=phi/*

/* Generate 3D lattice with 2D spins*/

#define N 16
static double lattice[N][N][N];
static double beta=0.454165;
static int nnx[N];
static int nny[N];
static int nnz[N];

double calc_energy();
double * calc_M();
void update();
double E_site(int i,int j,int k, double d);

void main() {
  int i,j,k;
  FILE *f;
  double E;
  double* M;
  init_genrand64(time(NULL));
  for (i=0; i<N; i++){
    for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
	lattice[i][j][k]=2*M_PI*genrand64_real3();
      }
    }
  }
  /* Generate the nn-arrays for every site*/
  for (int x=0; x < N; x++) {
    nnx[x]=(x+1) % N;
    nny[x]=(x+1) % N;
    nnz[x]=(x+1) % N;
  }
  
  f = fopen("Energies.txt","w");
  for (int a = 0; a<1000; a++) {
    update();
    E=calc_energy();
    M=calc_M();
    fprintf(f,"%f\t",E);
    fprintf(f,"%f\t",M[0]);
    fprintf(f,"%f\t",M[1]);
    fprintf(f,"%f\n",M[2]);
  }
  fclose(f);
}

void update() {
  double delta;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
	delta = 2*M_PI*genrand64_real3();
	double now=E_site(i,j,k,0);
	double trial=E_site(i,j,k,delta);
	if (now > trial) {
	  lattice[i][j][k]=trial;
	}
	else {
	  double p=genrand64_real3();
	  if (p<exp(-beta*(trial-now))) lattice[i][j][j]=lattice[i][j][k]+delta;
	}
      }
    }
  }	
}

double calc_energy() {
  /* check boundaries, implement periodic boundary conditions*/
  double E=0;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
	  E-=cos(lattice[i][j][k]-lattice[nnx[i]][j][k]);
	  E-=cos(lattice[i][j][k]-lattice[i][nny[j]][k]);
	  E-=cos(lattice[i][j][k]-lattice[i][j][nnz[k]]);
      }
    }
  }
  return E;
}
/*Laita tämä toimimaan*/
double* calc_M() {
  static double M[3];
  M[0]=0;
  M[1]=0;
  M[2]=0;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
	M[0]+=sin(lattice[i][j][k]);
	M[1]+=cos(lattice[i][j][k]);
	M[2]+=M[0]+M[1];
      }
    }
  }
  return M;
}


double E_site(int i, int j, int k, double d) {
  double E_ijk=0;
  /* E_ijk=cos(lattice[i][j][k]-lattice[i][j][k+1])+cos(lattice[i][j][k]-lattice[i][j+1][k])+cos(lattice[i][j][k]-lattice[i+1][j][k])+cos(lattice[i][j][k]-lattice[i-1][j][k])+cos(lattice[i][j][k]-lattice[i][j-1][k])+cos(lattice[i][j][k]-lattice[i][j][k-1]);*/

  if(i==N-1) {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[0][j][k])+cos(lattice[i][j][k]+d-lattice[i-1][j][k]);
  }
  else if(i==0) {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[N-1][j][k])+cos(lattice[i][j][k]+d-lattice[i+1][j][k]);
  }
  else {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[i+1][j][k])+cos(lattice[i][j][k]+d-lattice[i-1][j][k]);
  }

  if(j==N-1) {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[i][0][k])+cos(lattice[i][j][k]+d-lattice[i][j-1][k]);
  }
  else if(j==0) {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[i][N-1][k])+cos(lattice[i][j][k]+d-lattice[i][j+1][k]);
  }
  else {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[i][j+1][k])+cos(lattice[i][j][k]+d-lattice[i][j-1][k]);
  }

  if(k==N-1) {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[i][j][0])+cos(lattice[i][j][k]+d-lattice[i][j][k-1]);
  }
  else if(k==0) {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[i][j][N-1])+cos(lattice[i][j][k]+d-lattice[i][j][k+1]);
  }
  else {
    E_ijk+=cos(lattice[i][j][k]+d-lattice[i][j][k+1])+cos(lattice[i][j][k]+d-lattice[i][j][k-1]);
  }
  return -E_ijk;
}

