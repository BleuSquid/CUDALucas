#ifndef DRCSrw_h
# define DRCSrw_h "$Id: rw.h,v 8.1 2007/06/23 22:33:35 wedgingt Exp $"

# include <ctype.h>

# ifdef _MSC_VER
#  include <io.h>
# else
#  include <unistd.h>
# endif

# include "zero.h"

extern UL chkpnt_iterations;

extern int device_number; //msf

ATTRIBUTE_CONST const char *archive_name (void);
extern FILE *open_archive (void);
extern void print_time (int iterations = 0, int current = 0, int total = 0);

#   define FACTOR char
#   define EXPONENT UL
extern int check_point (UL q, UL n, UL j, double err, double *x);
extern void printbits(double *x, UL q, UL n, UL totalbits, UL b, UL c, double hi, double lo,
    const char *version_info, FILE *outfp, FILE *dupfp, int iterations, int current_iteration, bool archive = false);
extern int input (int argc, char **argv, EXPONENT *q, UL *n, UL *j, double *err, double **x,
                 EXPONENT last, char *M, FILE **infp, FILE **outfp, FILE **dupfp);
#endif /* ifndef DRCSrw_h */
