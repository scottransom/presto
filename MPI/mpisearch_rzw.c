#include "fftapps.h"
#include "mpi.h"
#ifdef LOGGING
#include "mpe.h"
#endif

typedef struct FFTINFO {
    double N;       /* Number of points in original time series      */
    double dt;      /* Sample duration per point of time series      */
    long nbins;     /* Number of points in the FFT file              */
    float nph;      /* Frequency 0 bin amplitude of the FFT          */
} fftinfo;

/* Some global datatypes */

MPI_Datatype fftinfotype;
MPI_Datatype positiontype;
MPI_Datatype fourierpropstype;

/* Function definitions */

void usage(void);
void make_position_struct(void);
void make_fftinfo_struct(void);
extern void master(int numprocs, int argc, char *argv[]);
extern void slave(int myid, int numprocs, int argc, char *vargv[]);

int main(int argc, char *argv[])
{
  int myid, numprocs;
  double startwtime, endwtime;

  /* Initialize and determine our ranks */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
#ifdef LOGGING
  MPE_Init_log();
#endif
  
  /* Start everyone off synchronized */

  MPI_Barrier(MPI_COMM_WORLD);
  
#ifdef LOGGING
  MPE_Start_log();
#endif
  
  if (myid == 0) {

    /* show the usage note if needed */

    if ((argc < 3) || (argc > 7)) usage();

    /* the master process */

    startwtime = MPI_Wtime();
    master(numprocs, argc, argv);
    endwtime = MPI_Wtime();

    /* Output the elapsed wall time */

    printf("/nWall clock time = %f\n", endwtime-startwtime);         
      
  } else {
    
    /* the slave processes */

    slave(myid, numprocs, argc, argv);

  }
  
#ifdef LOGGING
  MPE_Finish_log("mpisearch_rzw_log");
#endif
  
  MPI_Finalize();
  return(0);
}
  
void make_fftinfo_struct(void)
{
  int blockcounts[3] = {2,1,1};
  MPI_Datatype types[3] = {MPI_DOUBLE, MPI_LONG, MPI_FLOAT};
  MPI_Aint displs[3];
  int i;
  fftinfo finf;

  MPI_Address(&finf.N, &displs[0]);
  MPI_Address(&finf.nbins, &displs[1]);
  MPI_Address(&finf.nph, &displs[2]);
  for (i = 0 ; i < 3 ; i++)
    displs[i] -= displs[0];
  MPI_Type_struct(3, blockcounts, displs, types, &fftinfotype);
  MPI_Type_commit(&fftinfotype);
}


void make_position_struct(void)
{
  int blockcounts[3] = {1,3,1};
  MPI_Datatype types[3] = {MPI_FLOAT, MPI_DOUBLE, MPI_UB};
  MPI_Aint displs[3];
  int i;
  position pos[2];

  MPI_Address(&pos[0].pow, &displs[0]);
  MPI_Address(&pos[0].p1, &displs[1]);
  MPI_Address(&pos[1].pow, &displs[2]);
  for (i = 0 ; i < 3 ; i++)
    displs[i] -= displs[0];
  MPI_Type_struct(3, blockcounts, displs, types, &positiontype);
  MPI_Type_commit(&positiontype);
}


void make_fourierpropstype_struct(void)
{
  int blockcounts[7] = {1,1,1,1,1,12,1};
  MPI_Datatype types[7] = {MPI_DOUBLE, MPI_FLOAT, MPI_DOUBLE, MPI_FLOAT, \
			   MPI_DOUBLE, MPI_FLOAT, MPI_UB};
  MPI_Aint displs[7];
  int i;
  fourierprops props[2];

  MPI_Address(&props[0].r, &displs[0]);
  MPI_Address(&props[0].rerr, &displs[1]);
  MPI_Address(&props[0].z, &displs[2]);
  MPI_Address(&props[0].zerr, &displs[3]);
  MPI_Address(&props[0].w, &displs[4]);
  MPI_Address(&props[0].werr, &displs[5]);
  MPI_Address(&props[1].r, &displs[6]);
  for (i = 0 ; i < 7 ; i++)
    displs[i] -= displs[0];
  MPI_Type_struct(7, blockcounts, displs, types, &fourierpropstype);
  MPI_Type_commit(&fourierpropstype);
}


void usage(void)
{
    printf("\nUsage:  mpisearch_rzw file ncand [zlo] [lofreq] [rlo] [rhi]\n\n");
    printf("  Mandatory arguments:\n");
    printf("       'file' = a string containing the FFT file's name.\n");
    printf("                You must have an '.inf' file of the same\n");
    printf("                name as well.  Do not add a suffix.\n");
    printf("                Candidates will be returned in a file called\n");
    printf("                'file_rzw.cand'.  A Postscript format\n");
    printf("                candidate list will be in 'file_rzw.ps'.\n");
    printf("      'ncand' = (int) The routine will return 'ncand' candidates.\n");
    printf("                Must be less than or equal to 5000.\n");
    printf("  Optional arguments:\n");
    printf("        'zlo' = (int) Lowest Fourier freq deriv to search.\n");
    printf("     'lofreq' = (int) Lowest Fourier bin in FFT file.  This is\n");
    printf("                useful if you chop a long FFT for space reasons.\n");
    printf("                If 'lofreq' is present and not equal to '0', \n");
    printf("                we will assume the 0th frequency bin = 1 for\n");
    printf("                normalization purposes.\n");
    printf("        'rlo' = (int) lowest Fourier bin to search.\n");
    printf("        'rhi' = (int) highest Fourier bin to search.\n\n");
    printf("  'mpisearch_rzw' will search a region of the f-fdot plane for \n");
    printf("  pulsations in a file containing a long, single precision FFT\n");
    printf("  using the Correlation method (i.e. Ransom and \n");
    printf("  Eikenberry, 1997, unpublished as of yet).\n");
    printf("  The search uses a spacing of 0.5 frequency bins in\n");
    printf("  the fourier frequency (r) direction and 2 'bins' in\n");
    printf("  the fdot (z) direction.\n");
    printf("  The routine outputs formatted statistics for the 'ncand'\n");
    printf("  best candidates found in the search region.  If the\n");
    printf("  optional arguments are ommitted, the routine will search\n");
    printf("  the whole FFT file and assume the first bin is freq=0.\n\n");
    printf("  The routine currently cannot search the 'w' dimension.\n");
    printf("  This may be fixed shortly.\n");
    printf("                                        7 Nov 1997\n\n");
}

