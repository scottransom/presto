#define PREPGEOM   "mkdir /tmp/geomview"
#define PLOT3DFILE "/tmp/geomview/OOGL"
#include <stdio.h>

typedef struct COLOR {
    float r, g, b;
} color;

void threedplot(float **array, int nx, int ny, char *palette);
void plot3d_complex(float **array, int nx, int ny, char *palette);
void quick3d_complex(FILE * fp, double fftfreq, char *palette);
void hue(color *x, double val);
void heat(color *x, double val);
void rainbow(color *x, double val);
void gray(color *x, double val);
