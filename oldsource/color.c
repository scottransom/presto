#include "plot3d.h"

void hue(color *x, double val)
{
  double h, v = 1.0, r = 0, g = 0, b = 0, f, p, q, t;
  int i;

  if (val == 1.0)
    val = 0.0;
  h = val * 6.0;
  i = ((int) h > h) ? (int) h - 1 : (int) h;
  f = h - i;
  p = 0.0;
  q = 1 - f;
  t = f;
  switch (i) {
  case 0:
    r = v;
    g = t;
    b = p;
    break;
  case 1:
    r = q;
    g = v;
    b = p;
    break;
  case 2:
    r = p;
    g = v;
    b = t;
    break;
  case 3:
    r = p;
    g = q;
    b = v;
    break;
  case 4:
    r = t;
    g = p;
    b = v;
    break;
  case 5:
    r = v;
    g = p;
    b = q;
    break;
  }
  x->r = r;
  x->g = g;
  x->b = b;
}


void heat(color *x, double val)
{
  /*  red  */
  if (val < 0.70) {
    x->r = 1.42857 * val;
  } else
    x->r = 1.0;

  /*  green  */
  if (val < 0.48) {
    x->g = 0.0;
  } else
    x->g = 1.92307 * val - 0.92307;

  /*  blue  */
  if (val < 0.75) {
    x->b = 0.0;
  } else
    x->b = 4.0 * val - 3.0;
}


void rainbow(color *x, double val)
{
  /*  red  */
  if (val < 0.225) {
    x->r = 1.0;
  } else if (val < 0.4) {
    x->r = -5.71428 * val + 2.285712;
  } else if (val < 0.775) {
    x->r = 0.0;
  } else if (val < 0.965) {
    x->r = 5.26315 * val - 4.078941;
  } else
    x->r = 1.0;

  /*  green  */
  if (val < 0.015) {
    x->g = 0.0;
  } else if (val < 0.225) {
    x->g = 4.76190 * val - 0.0714285;
  } else if (val < 0.59) {
    x->g = 1.0;
  } else if (val < 0.775) {
    x->g = -5.40540 * val + 4.189185;
  } else
    x->g = 4.44444 * val - 3.44444;

  /*  blue  */
  if (val < 0.4) {
    x->b = 0.0;
  } else if (val < 0.6) {
    x->b = 5.0 * val - 2.0;
  } else if (val < 0.955) {
    x->b = 1.0;
  } else
    x->b = -4.861111 * val + 5.642361;
}


void gray(color *x, double val)
{
  /*  red, green, and blue  */
  x->r = x->g = x->b = val;
}
