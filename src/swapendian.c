#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tmpswap=(a);(a)=(b);(b)=tmpswap;
#endif

static unsigned char tmpswap;

float swap_float(float var)
{
    unsigned char *buffer;
    float *fptr;

    buffer = (unsigned char *) (&var);
    SWAP(buffer[0], buffer[3]);
    SWAP(buffer[1], buffer[2]);
    fptr = (float *) buffer;
    return *fptr;
}

double swap_double(double var)
{
    unsigned char *buffer;
    double *dptr;
    buffer = (unsigned char *) (&var);
    SWAP(buffer[0], buffer[7]);
    SWAP(buffer[1], buffer[6]);
    SWAP(buffer[2], buffer[5]);
    SWAP(buffer[3], buffer[4]);
    dptr = (double *) buffer;
    return *dptr;
}

long double swap_longdouble(long double var)
{
    unsigned char *buffer;
    long double *ldptr;

    buffer = (unsigned char *) (&var);
    SWAP(buffer[0], buffer[11]);
    SWAP(buffer[1], buffer[10]);
    SWAP(buffer[2], buffer[9]);
    SWAP(buffer[3], buffer[8]);
    SWAP(buffer[4], buffer[7]);
    SWAP(buffer[5], buffer[6]);
    ldptr = (long double *) buffer;
    return *ldptr;
}

long long swap_longlong(long long var)
{
    unsigned char *buffer;
    long long *llptr;

    buffer = (unsigned char *) (&var);
    SWAP(buffer[0], buffer[7]);
    SWAP(buffer[1], buffer[6]);
    SWAP(buffer[2], buffer[5]);
    SWAP(buffer[3], buffer[4]);
    llptr = (long long *) buffer;
    return *llptr;
}

int swap_int(int var)
{
    unsigned char *buffer;
    int *iptr;

    buffer = (unsigned char *) (&var);
    SWAP(buffer[0], buffer[3]);
    SWAP(buffer[1], buffer[2]);
    iptr = (int *) buffer;
    return *iptr;
}

unsigned int swap_uint(unsigned int var)
{
    unsigned char *buffer;
    unsigned int *iptr;

    buffer = (unsigned char *) (&var);
    SWAP(buffer[0], buffer[3]);
    SWAP(buffer[1], buffer[2]);
    iptr = (unsigned int *) buffer;
    return *iptr;
}

short swap_short(short var)
{
    unsigned char *buffer;
    short *sptr;

    buffer = (unsigned char *) (&var);
    SWAP(buffer[0], buffer[1]);
    sptr = (short *) buffer;
    return *sptr;
}

unsigned short swap_ushort(unsigned short var)
{
    unsigned char *buffer;
    unsigned short *sptr;

    buffer = (unsigned char *) (&var);
    SWAP(buffer[0], buffer[1]);
    sptr = (unsigned short *) buffer;
    return *sptr;
}

#undef SWAP
