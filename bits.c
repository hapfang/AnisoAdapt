#include "bits.h"

int setbit (int n, int i)   /* to set the rightmost bit, use i=0 etc. */
{
  return n | (01 << i);
}

int getbit (int n, int i)
{
  return 01 & (n >> i);
}
