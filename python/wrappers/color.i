%module color

%{
#include <stdio.h>
#include "plot3d.h"
%}

typedef struct COLOR {
    float r, g, b;
} color;

void hue(color *x, float val);
void heat(color *x, float val);
void rainbow(color *x, float val);
void gray(color *x, float val);







