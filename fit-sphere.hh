#pragma once

// data : sample data points (x0, y0, z0, x1, y1, z1, ..., xm, ym, zm)
// m : number of data points
// R : fixed radius
// result : the approximated center point
// tolerance : end the iteration when the norm of the change is less than this
// max_iteration : maximum iteration count
// lambda : step length for the Newton method
// Returns true when the tolerance was achieved before the maximum iteration count
//   has been reached; false otherwise.
bool fit_sphere(const double *data, size_t m, double R, double *result,
                double tolerance = 1.0e-7, double max_iteration = 1000, double lambda = 0.1);
