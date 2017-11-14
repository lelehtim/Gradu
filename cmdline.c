/* Some command line processing tools */
/* void cmdl_init( int argc, char *argv[] )    
 *    inits tool
 *
 * int cmdl_get_opt()                 
 *    if next is option (-x), returns 'x', otherwise 0
 *
 * int cmdl_get_int( char * error )   
 *    next MUST be int, otherwise prints error and exits
 *
 * double cmdl_get_double( char *error )    
 *    same with double
 *
 * int cmdl_is_char( char c )
 *    if next character is c, returns 1 and skips c, otherwise returns 0
 *
 * char * cmdl_get_string( char *error ) 
 *    retuns next argument as a string, if no args prints error and exits
 *
 * int cmdl_is_args()
 *    returns the number of argv-elements remaining
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_ARG 200

void cmdl_init( int argc, char *argv[] );
int cmdl_get_opt();
int cmdl_get_int(char *err);
double cmdl_get_double(char *err);
int cmdl_is_char(char c);
char * cmdl_get_string( char *err );
int cmdl_args_remain();


static char *args[MAX_ARG];
static char **argp;  /* current arg-element */
static char *ap;     /* current char position in arg */
static int argcount, narg;

/* error handling */
void cmdl_err(char *err)
{
  fprintf(stderr," Error in command line arguments\n");
  fprintf(stderr,err);
  exit(0);
}


/* set up the arg variables - first is always the prog name! */

void cmdl_init( int argc, char *argv[] )
{
  int i;
  for (i=0; i<argc && i<MAX_ARG; i++) args[i] = argv[i];
  argcount = argc;

  if (argcount > 1) {
    argp = args+1;
    narg = argcount-1;
    ap   = *argp;
  } else {
    ap = NULL;  /* Null means nothing to do */
  }
}  

/* work function - get next item */

void cmdl_do_next_item()
{
  if (narg > 0) {
    argp++;
    narg--;
    ap = *argp;
  } else ap = NULL;
}

/* work function */

char * cmdl_find_next_string()
{
  if (ap == NULL) return(NULL);
  if (*ap == 0) cmdl_do_next_item();
  if (ap == NULL) return(NULL);
  return( ap );
}

/* return possible option char */

int cmdl_get_opt()
{
  int c;
 
  if (cmdl_find_next_string() == NULL) return(0);
  if (*ap != '-') return(0); /* not a - */
  ap++;
  if (*ap == 0) return('-'); /* just a '-' char */
  c = *ap;
  ap++;
  return( c );
}


/* command line processing: return int, or quit with error */
int cmdl_get_int(char *err)
{
  int t;
  char *ep;

  if (cmdl_find_next_string() == NULL) cmdl_err(err);
  
  t = strtol( ap, &ep, 10);
  if (ep == ap) cmdl_err(err);

  ap = ep; /* set pointer past the number */
  
  return( t );
}


/* command line processing: return double or exit */
double cmdl_get_double(char *err)
{
  double t;
  char *ep;

  if (cmdl_find_next_string() == NULL) cmdl_err(err);
  
  t = strtod( ap, &ep );
  if (ep == ap) cmdl_err(err);

  ap = ep; /* set pointer past the number */
  
  return( t );
}


/* command line processing: is next char c? if so, skip it */
int cmdl_is_char(char c)
{
  if (cmdl_find_next_string() == NULL || *ap != c) return(0);
  ap++;
  return(1);
}

/*  retuns next argument as a string, 
 *  if no args prints error and exits
 */ 

char * cmdl_get_string( char *err ) 
{
  char *p;

  if (cmdl_find_next_string() == NULL) cmdl_err(err);
  p = ap;
  cmdl_do_next_item();   /* go to next, so that won't return 
			  * this again 
			  */
  return( p );
}

/* return number of args remaining */

int cmdl_args_remain()
{
  if (cmdl_find_next_string() == NULL) return(0);

  /* now ap points to an argument string, and narg 
     counts the others - so return */
  return( narg + 1 );
}
