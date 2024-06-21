using System;

class NastyPSO {    
    static int dimension = 30;
    static int size = 20;
    static double a = 1.49180;
    static double w = 0.729844;
    static double iterations = 10000;
    static double repeats = 100;
    
    static double [][] position;
    static double [][] bestpos;
    static double [][] velocity;
    static double [] fitness;
    static double [] bestfit;
    static int gbest;
    static double range = 1000;
    
    static int idnum;
    static int inext;
    static int inextp;
    static long [] ma = new long[56];
    static int iff = 0;
    
    static long MBIG = 1000000000;
    static long MSEED = 161803398;
    static int MZ = 0;
    static double FAC = (1.0/MBIG);

    static double ran3() {
	long mj, mk;
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
    
    static double f(double [] x) {
	double tmp = 0.0;
	for (int i = 0; i < dimension; ++i) {
	    tmp += (x[i] - (i + 1)) * (x[i] - (i + 1));
	}
	return tmp;
    }
    
    static void init() {
	for (int i = 0; i < size; ++i) {
	    for (int j = 0; j < dimension; ++j) {
		position[i][j] = range * ran3() - range / 2;
		bestpos[i][j] = position[i][j];
		velocity[i][j] = 0;
	    }
	    fitness[i] = 1e300;
	    bestfit[i] = 1e300;
	}
	gbest = 0;
    }
    
    static void update(int k) {
	double v;
	for (int i = 0; i < dimension; ++i) {
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
	for (int i = 0; i < dimension; ++i) {
	    position[k][i] += velocity[k][i];
	}
    }
    
    public static void Main() {
	Random r = new Random();
	idnum = - r.Next(2147483647);
	position = new double [size][];
	bestpos = new double [size][];
	velocity = new double [size][];
	for (int i = 0; i < size; ++i) {
	    position[i] = new double[dimension];
	    bestpos[i] = new double[dimension];
	    velocity[i] = new double[dimension];
	}
	fitness = new double[size];
	bestfit = new double[size];
	
	for (int i = 0; i < repeats; ++i) {
	    init();
	    for (int j = 0; j < iterations; ++j) {
		for (int k = 0; k < size; ++k) {
		    fitness[k] = f(position[k]);
		    if (fitness[k] < bestfit[k]) {
			bestfit[k] = fitness[k];
			for (int l = 0; l < dimension; ++l) {
			    bestpos[k][l] = position[k][l];
			}
			if (bestfit[k] < bestfit[gbest]) {
			    gbest = k;
			}
		    } 
		    update(k);
		}
	    }
	    Console.Write(i + " : " + bestfit[gbest] + " [");
	    for (int j = 0; j < dimension - 1; ++j) {
		Console.Write(position[gbest][j] + ", ");
	    }
	    Console.WriteLine(position[gbest][dimension - 1] + "]");
	}
    }   
};
