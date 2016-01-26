#include <stdlib.h>
#include <math.h>

/* max(A,B) - larger (most +ve) of two numbers (generic) */
#define max(A,B) ((A)>(B)?(A):(B))
/* min(A,B) - smaller (least +ve) of two numbers (generic) */
#define min(A,B) ((A)<(B)?(A):(B))


void apprgrdn(unsigned short n,
              double g[],
              double x[],
              double f, double fun(), double deltax[], unsigned short obj)
{
/* Function APPRGRDN performs the finite difference approximation 
   of the gradient <g> at a point <x>.
   f      is the calculated function value at a point <x>,
   <fun>  is the name of a function that calculates function values,
   deltax is an array of the relative stepsizes.
   obj    is the flag indicating whether the gradient of the objective
          function (1) or the constraint function (0) is to be calculated. 
*/
    double const lowbndobj = 2.0e-10, lowbndcnt = 5.0e-15, ten = 10.0, half = 0.5;
    double d, y, fi;
    unsigned short i, j, center = 0;
    for (i = 0; i < n; i++) {
        y = x[i];
        d = max(lowbndcnt, fabs(y));
        d *= deltax[i];
        if (obj) {
            if (fabs(d) < lowbndobj) {
                if (deltax[i] < 0.0)
                    d = -lowbndobj;
                else
                    d = lowbndobj;
                center = 1;
            } else
                center = 0;
        } else if (fabs(d) < lowbndcnt) {
            if (deltax[i] < 0.0)
                d = -lowbndcnt;
            else
                d = lowbndcnt;
        }
        x[i] = y + d;
        fi = fun(x);
        if (obj) {
            if (fi == f) {
                for (j = 1; j <= 3; j++) {
                    d *= ten;
                    x[i] = y + d;
                    fi = fun(x);
                    if (fi != f)
                        break;
                }
            }
        }
        g[i] = (fi - f) / d;
        if (obj) {
            if (center) {
                x[i] = y - d;
                fi = fun(x);
                g[i] = half * (g[i] + (f - fi) / d);
            }
        }
        x[i] = y;
    }
}
