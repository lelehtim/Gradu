/*************************************************************************
 *
 * Analyses data for errors and autocorrelations
 *
 * COMPILE:
 *  cc -o errors errors.c cmdline.c -lm
 *
 * Use:
 * errors [-opt] file

 *  Reads in file and prints out average, error and autocorrelation time
 *  from the "file".  This is assumed to have measurements so that
 *  each line corresponds to measurements on one configuration.  Different
 *  quantities are always on different columns (separated by one or more
 *  spaces or tab-characters).  Program automatically
 *  recognizes rows and columns in the file.  Example:
 *  
 *  a_1 b_1 c_1 
 *  a_2 b_2 c_2 
 *  a_3 b_3 c_3 
 *  ...
 *
 * Options:
 * -c c1,c2,...  : Comma-separated list of columns selected for
 *                 analysis. First column is 1
 *                 At least one column must be selected!
 * -n n          : Use (at most) n measurements in the analysis
 * -s n          : Skip n measurements from beginning
 * -b n          : Error analysis by blocking the  measurements to 
 *                 block length n.
 *                 This is alternative to autocorrelation analysis, so
 *                 no autocorrelations printed now.
 * -t            : do not calculate autocorrelations, i.e. use only
 *                 naive errors
 * -T length     : Instead of error analysis, print out the autocorrelation
 *                 function up to lag length.
 *
 * Examples:
 *
 * > errors -c3 meas 
 * Prints out average, error and autocorrelation time of the data
 * in column 3.
 *
 * > errors -c3,4  -s2000 meas    
 * prints out average, error and autocorrelation time of cols 3 and 4,
 * after skipping first 2000 iterations.
 *
 * > errors -c3,4 -b500 -s2000 -i20000  meas    
 * Do now error analysis using blocks of 500 measurements, after
 * skipping 2000 measurements from the beginning and including only
 * 20000 following measurements.
 *
 *************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



/* Prototypes */
int readdata(FILE *f, double *d[], int cols[], int ncol, int skip, int nmeas);

void errorcalc(double d[], int n, int autocorr,
	       double *avep, double *sigp, double *tintp);
void blockerr(double d[], int n, int blocksize, double *avep, double *sigp);
void autocorrfunc(double d[], int n, int maxlag);

/* command line handling - see cmdline.c */
void cmdl_init( int argc, char *argv[] );
int cmdl_get_opt();
int cmdl_get_int(char *err);
double cmdl_get_double(char *err);
int cmdl_is_char(char c);
char * cmdl_get_string( char *err );
int cmdl_args_remain();



#define sqr(a) ((a)*(a))

char usage[]  = " Use: errors [-opt] file\n\
 Reads in file and prints out average, error and autocorrelation time\n\
  c c1,c2,...  : comma-separated list of columns, first is 1\n\
                 at least one column is required\n\
  n n          : use only n measurements\n\
  s n          : skip n measurements from beginning\n\
  b n          : block measurements to block length n, no autocorrelations\n\
  t            : no autocorrelations\n\
  T length     : print autocorrelation function up to distance length\n";


/* This routine prints usage string and exits */
void err_args()
{  
  fprintf(stderr,usage);
  exit(0); /* quit the program */
}


#define MAXCOLS 30

int main(int argc,char * argv[])
{
  int col[MAXCOLS]; /* data columns */
  int ncols = 0;    /* number of columns */
  int nmeas = 0;    /* use all meas */
  int skip  = 0;    /* do not skip */
  int autoc = 1;    /* measure autocorrelations */
  int autof = 0;    /* do not print autocorr. function */
  int block = 0;    /* no blocking */
  FILE *f;          /* datafile */
  double *d[MAXCOLS]; /* data arrays */
  int n;            /* data array size */
  int i,opt;
  char *fname;      /* data file name */

  /* Process all command line arguments and filename */

  if (argc <= 1) err_args();

  /* command line arguments (-) */
  cmdl_init( argc, argv );

  while (opt = cmdl_get_opt()) {
    switch(opt) {

    case 't': autoc = 0; break;
    case 'T': autof = cmdl_get_int(usage); break;
    case 'n': nmeas = cmdl_get_int(usage); break;
    case 's': skip  = cmdl_get_int(usage); break;
    case 'b': block = cmdl_get_int(usage); break;
    case 'c': 
      do {
	if (ncols >= MAXCOLS) {
	  fprintf(stderr,"Too many columns (max %d)\n",MAXCOLS);
	  exit(0);
	}
	col[ncols++] = cmdl_get_int(usage);
      } while(cmdl_is_char(','));
      break;

    default: err_args();
    } /* switch */
  } /* while */

  if (ncols == 0) err_args();    /* col has to be given */
  fname = cmdl_get_string(usage);    /* get filename */

  /* Now *argv points to the file name */
  f = fopen(fname,"r");
  if (f == NULL) {
    fprintf(stderr," Could not open file %s\n",fname);
    exit(0);
  }
  
  /* read in the data - n: number of elements read */
  n = readdata( f, d, col, ncols, skip, nmeas );

  for (i=0; i<ncols; i++) {
    if (autof) autocorrfunc(d[i],n,autof);
    else {
      double ave,sig,tint;

      if (block) {
	blockerr(d[i],n,block,&ave,&sig);
	
	if (i==0) printf(" col - average - sigma\n");
	printf(" %d  %.10g  %.8g\n",col[i],ave,sig);
      
      } else {
	errorcalc(d[i],n,autoc,&ave,&sig,&tint);
	
	if (i==0) printf(" col - average - sigma - tau_int\n");
	printf(" %d  %.10g  %.8g  %g\n",col[i],ave,sig,tint);
      }
    }
  }

  return(1);
}


/* calculate and print average, errors and autocorrelation time 
 * Calculate integrated 
 *   tint = 0.5 + sum_{t=1}^N C(t)
 * where N-sum is stopped when N >= TINTSTOP*tint
 *   C(t) = < (d[i]-<d>) (d[i-t] - <d>)> / <(d-<d>)^2>
 */
#define TINTSTOP 6

void errorcalc(double d[], int n, int autocorr,
	       double *avep, double *sigp, double *tintp)
{
  double ave,sig,tint;
  int i;

  ave = sig = 0.0;

  for(i=0; i<n; i++) ave += d[i];
  ave /= n;
  for(i=0; i<n; i++) sig += sqr(d[i]-ave);
  sig /= n;

  /* time correlations */

  if (autocorr) {
    int it;     /* autocorrelation lag */

    tint = 0.5;    
    /* start loop over lags */
    for (it=1; it < TINTSTOP*tint && it < n/2; it++) {
      double av1=0,av2=0;  /* local averages */
      double fi=0;   
      int nc = n-it;       /* number of measurements */
      int j;

      for (j=0; j<nc; j++) {
	fi  += d[j]*d[j+it];
	av1 += d[j];
	av2 += d[j+it];
      }

      fi = ( fi/nc - (av1/nc)*(av2/nc) )/( sig );
      tint += fi*nc/n;
    } 
    if (it >= n/2) fprintf(stderr," ** correlation > N/2*%d\n",TINTSTOP);
  }

  *avep = ave;
  if (autocorr) {
    *sigp  = sqrt( 2*fabs(tint) * sig/(n-1) );
    *tintp = tint;
  } else {
    *sigp  = sqrt( sig/(n-1) );
  }
}

/* Calculate errors with blocking
 * blocks data to blocksize-sized blocks; with
 * n/blocksize blocks; if division is not even rest not used.
 */

void blockerr(double d[], int n, int blocksize, double *avep, double *sigp)
{
  double ave,sig,*blockave;
  int i,id,iblock,nblocks;

  ave = sig = 0.0;
  nblocks = n/blocksize;  
  blockave = (double *)malloc(nblocks * sizeof(double));

  id = 0;
  for (iblock=0; iblock<nblocks; iblock++) {
    blockave[iblock] = 0;
    for(i=0; i<blocksize; i++,id++) blockave[iblock] += d[id];
    blockave[iblock] /= blocksize;
  }

  /* sum ave and sigma, assuming blockave is measurement */
  for (iblock=0; iblock<nblocks; iblock++) 
    ave += blockave[iblock];
  ave /= nblocks;
  for (iblock=0; iblock<nblocks; iblock++) 
    sig += sqr( blockave[iblock] - ave );
  sig = sqrt( sig/(nblocks*(nblocks-1)) );

  free( blockave );
  *avep = ave;
  *sigp = sig;
}

/* Print autocorrelation function to lat maxlag
 */

void autocorrfunc(double d[], int n, int maxlag)
{
  int it;     /* autocorrelation lag */
  int i;
  double ave=0,sig=0;

  for(i=0; i<n; i++) ave += d[i];
  ave /= n;
  for(i=0; i<n; i++) sig += sqr(d[i]-ave);
  sig /= n;

  fprintf(stderr,"lag -- autocorrelation\n");

  /* loop over lags */
  for (it=0; it <= maxlag && it < n/2; it++) {
    double av1=0,av2=0;  /* local averages */
    double fi=0;   
    int nc = n-it;       /* number of measurements */
    int j;

    for (j=0; j<nc; j++) {
      fi  += d[j]*d[j+it];
      av1 += d[j];
      av2 += d[j+it];
    }

    fi = ( fi/nc - (av1/nc)*(av2/nc) )/( sig );

    printf("%d %g\n",it,fi);

  } 
}  


/* readdata() 
 * reads in data from FILE *f, from columns given in array col[].
 * It skips skip lines from beginning
 * If nmax > 0 reads at most nmax lines (after skipping)
 * Returns the number of measurements, and
 * in d the data arrays.  It allocates the arrays, so pass
 * in only array of pointers (or address of a pointer)!
 * double *d[] must have ncols elements!
 *
 *   double *d;
 *   int col = 1;
 *   n = readdata(f, &d, &col, 1, 0, 0)  
 * reads col 1 to d[], n elements
 *
 *   double *d[2];
 *   int cols[] = {2,3};
 *   n = readdata(f, d, cols, 2, 0, 0)   
 *  reads cols 2 and 3 into d[0][i] and d[1][i], i=0..(n-1),
 *  respectively
 */
#define BUFSIZE 2000

int readdata( FILE *f,     /* input file */
	      double *d[], /* data arrays to be read */
	      int col[], int ncols,  /* columns to be read */
	      int skip,    /* how many to skip from beginning */
	      int nmax)    /* max number of measurements to be read */
{
  int datsize;      /* array size */
  char buf[BUFSIZE];   /* line buffer */
  int line;         /* line number */
  double *cdata;    /* data column buffer */
  int maxcol;       /* max column needed */
  int i,j;          /* work */

  /* check how many lines file has - if nmax, read only necessary lines */
  for (i=0; (!nmax || i < nmax+skip); i++)
    if (fgets( buf, BUFSIZE, f ) == NULL) break;

  datsize = i - skip;
  if (datsize <= 0) {
    fprintf(stderr," No data to be read\n");
    exit(0);
  }
  
  fprintf(stderr,"  Reading in %d columns, %d measurements\n",
	  ncols,datsize);

  /* allocate arrays */
  for (i=0; i<ncols; i++)
    d[i] = (double *)malloc( datsize * sizeof(double) );

  /* what is the max column we need, allocate work */
  for (maxcol=i=0; i<ncols; i++) {
    if (maxcol < col[i]) maxcol = col[i];
    if (col[i] <= 0) {
      fprintf(stderr," Not valid column %d\n",col[i]);
      exit(0);
    }
  }
  cdata = (double *)malloc( maxcol * sizeof(double) );

  /* back to beginning of file, start processing it */
  rewind(f);

  /* skip lines */
  for (line=0; line<skip; line++ ) fgets( buf, BUFSIZE, f );
  
  /* read in lines */
  for (i=0; i<datsize; i++) {
    char *tp,*cp; 

    line++;  /* line counter */
    fgets( buf, BUFSIZE, f );

    /* process the line - get columns using strtod() */
    cp = buf;
    for (j=0; j<maxcol; j++) {
      cdata[j] = strtod(cp, &tp);
      if (tp == cp) {
	fprintf(stderr," Line %d: not %d columns\n",line,maxcol);
	exit(0);
      }
      cp = tp;
    }
    /* now cdata contains the numbers, copy to d */
    for (j=0; j<ncols; j++) 
      d[j][i] = cdata[col[j]-1];
  } 
  
  /* all read, free cdata */
  free( cdata );
  return( datsize );
}

