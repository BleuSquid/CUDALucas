#include "setup.h"
#include "balance.h"
#include "globals.h"

UL is_big(UL j, UL big, UL small, UL n)
{
  return((((big*j) % n) >= small) || j == 0);
}

void balancedtostdrep(double *x, UL n, UL b, UL c, double hi, double lo, UL mask, UL shift)
{
  UL sudden_death = 0, j = 0, NminusOne = n - 1, k, k1;

  while(true)
  {
    k = j + ((j & mask) >> shift);
    if (x[k] < 0.0)
    {
      k1 = (j + 1) % n;
      k1 += (k1 & mask) >> shift;
      --x[k1];
      if (j == 0 || (j != NminusOne && is_big(j, b, c, n)))
        x[k] += hi;
      else
        x[k] += lo;
    }
    else
      if (sudden_death)
        break;
    if (++j == n)
    {
      sudden_death = 1;
      j = 0;
    }
  }
}
