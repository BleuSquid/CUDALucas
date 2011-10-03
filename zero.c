#include "setup.h"
#include "zero.h"

int is_zero(double *x, UL n, UL mask, UL shift)
{
  register UL j, offset;

  for (j = 0; j < n; ++j)
  {
    offset = j + ((j & mask) >> shift);
    if (rint(x[offset]))
      return(0);
  }
  return(1);
}
