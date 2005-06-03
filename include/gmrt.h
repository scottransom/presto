void set_GMRT_padvals(float *fpadvals, int good_padvals);
void get_GMRT_file_info(FILE *files[], char *datfilenms[], int numfiles, 
			float clipsig, long long *N, int *ptsperblock, int *numchan, 
			double *dt, double *T, int output);
void GMRT_update_infodata(int numfiles, infodata *idata);
int skip_to_GMRT_rec(FILE *infiles[], int numfiles, int rec);
int read_GMRT_rawblock(FILE *infiles[], int numfiles, unsigned char *data, 
		       int *padding);
int read_GMRT_rawblocks(FILE *infiles[], int numfiles, unsigned char rawdata[], 
			int numblocks, int *padding);
int read_GMRT(FILE *infiles[], int numfiles, float *data, int numpts, 
	      double *dispdelays, int *padding, int *maskchans, 
	      int *nummasked, mask *obsmask);
void get_GMRT_channel(int channum, float chandat[], 
		      unsigned char rawdata[], int numblocks);
int prep_GMRT_subbands(unsigned char *rawdata, float *data, double *dispdelays, 
		       int numsubbands, int transpose, int *maskchans, 
		       int *nummasked, mask *obsmask);
int read_GMRT_subbands(FILE *infiles[], int numfiles, float *data, 
		       double *dispdelays, int numsubbands, int transpose, 
		       int *padding, int *maskchans, int *nummasked, 
		       mask *obsmask);
void GMRT_hdr_to_inf(char *datfilenm, infodata *idata);
void convert_GMRT_block(short *indata, unsigned char *outdata);
void convert_GMRT_block_8bit(unsigned char *indata, unsigned char *outdata);
