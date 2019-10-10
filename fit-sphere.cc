#include <Eigen/Dense>

using namespace Eigen;

bool fit_sphere(const double *data, size_t m, double R, double *result,
                double tolerance, double max_iteration, double lambda) {
  const Map<const MatrixXd> X(data, 3, m);
  Map<Vector3d> a(result);
  a.Zero();
  for (size_t i = 0; i < m; ++i)
    a += X.col(i);
  a /= m;
  MatrixXd J(3 * m, 3);
  VectorXd rhs(3 * m);
  auto id = Matrix3d::Identity();
  for (size_t it = 0; it < max_iteration; ++it) {
    for (size_t i = 0; i < m; ++i) {
      auto v = X.col(i) - a;
      auto d = v.norm();
      J.block<3,3>(3*i,0) = id - (id - v * v.transpose() / (d * d)) * R / d;
      rhs.segment<3>(3*i) = v / d * (d - R);
    }
    Vector3d delta = J.householderQr().solve(rhs);
    a += lambda * delta;
    if (delta.norm() < tolerance)
      return true;
  }
  return false;
}
