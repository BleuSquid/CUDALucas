#define DRCSglobals_h "$Id: globals.h,v 8.1 2007/06/23 22:33:35 wedgingt Exp $"

extern double *cn, *sn, *two_to_phi, *two_to_minusphi, *scrambled;
extern double high, low, highinv, lowinv;
extern long *permute;
extern int partial; /* read a partial result? */
#ifdef DEBUG
extern double highest_err;
#endif
