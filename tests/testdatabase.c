#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../../include/presto.h"

int main(int argc, char *argv[])
{
  psrdata psr, *psrs;
  psrdatabase pdata;
  int np, pnum;

  np = read_database(&pdata);
  psrs = collect_psrdata(&np);
  psr = get_psrdata(argv[1]);
  printf("There are %d pulsars in the database.\n\n",np);
  printf(" RA = %f\n", psr.ra2000);
  printf("DEC = %f\n", psr.dec2000);
  printf("  p = %f\n", psr.p);
  pnum = psr_number_from_name(argv[1], &pdata);
  printf("\nPSR J%.12s is pulsar number %d\n\n", psr.jname, pnum);
  printf(" RA = %f\n", psrs[pnum].ra2000);
  printf("DEC = %f\n", psrs[pnum].dec2000);
  printf("  p = %f\n", psrs[pnum].p);
  free(psrs);
  return 1;
}
