#include <cstdlib>
#include <cstdio>
#include <ctime>

using namespace std;

unsigned int dimension = 30;
unsigned int size = 20;
double a = 1.49180;
double w = 0.729844;
double iterations = 10000;
double repeats = 100;

double **position;
double **bestpos;
double **velocity;
double *fitness;
double *bestfit;
unsigned int gbest;
double range = 1000;

int idnum;
int inext;
int inextp;
long ma[56];
int iff = 0;

#define MBIG   1000000000
#define MSEED  161803398
#define MZ     0
#define FAC    (1.0/MBIG)

double ran3(void) {
  register long mj, mk;
  int i, ii, k;

  if (idnum < 0 || iff == 0) {
    iff = 1;
    mj = MSEED - (idnum < 0 ? -idnum : idnum);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i=1; i <= 54; i++) {
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ) mk += MBIG;
      mj = ma[ii];
    }
    for (k=1; k <= 4; k++) {
      for (i=1; i <= 55; i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
      }
    }
    inext  = 0;
    inextp = 31;
    idnum  = 1;
  }
  // start of actual algorithm
  if (++inext == 56)  inext = 1;
  if (++inextp == 56) inextp = 1;
  mj = ma[inext] - ma[inextp];
  while (mj < MZ) mj += MBIG;
  while (mj > MBIG) mj -= MBIG;
  ma[inext] = mj;
  return mj*FAC;
}

double f(double *x) {
  double tmp = 0.0;
  for (unsigned int i = 0; i < dimension; ++i) {
    tmp += (x[i] - (i + 1)) * (x[i] - (i + 1));
  }
  return tmp;
}

void init() {
  for (unsigned int i = 0; i < size; ++i) {
    for (unsigned int j = 0; j < dimension; ++j) {
      position[i][j] = range * ran3() - range / 2;
      bestpos[i][j] = position[i][j];
      velocity[i][j] = 0;
    }
    fitness[i] = 1e300;
    bestfit[i] = 1e300;
  }
  gbest = 0;
}

void update(unsigned int k) {
  double v;
  for (unsigned int i = 0; i < dimension; ++i) {
    v = w * velocity[k][i] 
      + (bestpos[k][i] - position[k][i]) * a * ran3()
      + (bestpos[gbest][i] - position[k][i]) * a * ran3();
    if (v > range) {
      v = range;
    }
    else if (v < -range) {
      v = -range;
    }
    velocity[k][i] = v;
  }
  for (unsigned int i = 0; i < dimension; ++i) {
    position[k][i] += velocity[k][i];
  }
}

int main(int argc, char **argv) {
  srand(time(0));
  idnum = -rand();
  position = new double *[size];
  bestpos = new double *[size];
  velocity = new double *[size];
  for (unsigned int i = 0; i < size; ++i) {
    position[i] = new double[dimension];
    bestpos[i] = new double[dimension];
    velocity[i] = new double[dimension];
  }
  fitness = new double[size];
  bestfit = new double[size];

  for (unsigned int i = 0; i < repeats; ++i) {
    init();
    for (unsigned int j = 0; j < iterations; ++j) {
      for (unsigned int k = 0; k < size; ++k) {
	fitness[k] = f(position[k]);
	if (fitness[k] < bestfit[k]) {
	  bestfit[k] = fitness[k];
	  for (unsigned int l = 0; l < dimension; ++l) {
	    bestpos[k][l] = position[k][l];
	  }
	  if (bestfit[k] < bestfit[gbest]) {
	    gbest = k;
	  }
	} 
	update(k);
      }
    }
    printf("%d : %f [", i, bestfit[gbest]);
    for (unsigned int j = 0; j < dimension - 1; ++j) {
      printf("%f, ", position[gbest][j]);
    }
    printf("%f]\n", position[gbest][dimension - 1]);
  }

  for (unsigned int i = 0; i < size; ++i) {
    delete [] position[i];
    delete [] bestpos[i];
    delete [] velocity[i];
  }
  delete [] fitness;
  delete [] bestfit;

  return 0;
}
