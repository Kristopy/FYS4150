#include "sym_tridiag_solver.h"

void sym_tridiag_solver(vec &x, vec y, int N, double b, double d)
{

  /*

    | 1            |    |x0  |   |y0  |
    | b d b        |    |x1  |   |y1  |
    |     ...      | *  |... | = |... |
    |         b d b|    |    |   |    |
    |             1|    |xN-1|   |yN-1|

  */
  x.set_size(N);

  double bb = b*b;

  vec diag(N); diag.fill(d);

  y(1) += -b*y(0);
  for (int i = 2; i < N - 1; ++i) {
    diag(i) += -bb/diag(i-1);
       y(i) += -b*y(i-1)/diag(i-1);
  }

  x(0) = y(0); x(N - 1) = y(N - 1);
  for (int i = N - 2; i > 0; --i) {
    x(i) = (y(i) - b*x(i + 1))/diag(i);
  }

  return;
}
