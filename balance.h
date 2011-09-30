#define DRCSbalance_h "$Id: balance.h,v 8.1 2007/06/23 22:33:35 wedgingt Exp $"

extern void balancedtostdrep (double *x, UL n, UL b, UL c, double hi, double lo, UL mask, UL shift);
extern UL is_big (UL j, UL big, UL small, UL n);
