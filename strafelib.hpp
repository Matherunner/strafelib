/** @file */

#include <cmath>

/// The initial jumping speed, before gravity is applied.
constexpr const double JUMP_SPEED = 268.3281572999748;

/// Compute the speed after applying ground friction.
///
/// This function runs at constant time. The caller is responsible of computing
/// \p tau_k correctly by considering the values of the entity friction and
/// the edgefriction.
///
/// If you are working on squared speeds for performance reasons, use fric_speedsq()
/// instead to avoid computing square roots for \p speed when the geometric friction
/// is in effect.
double fric_speed(double speed, double E, double tau_k)
{
    if (speed >= E) {
        return speed * (1 - tau_k);
    }

    const double tau_E_k = tau_k * E;
    if (speed >= tau_E_k && speed >= 0.1) {
        return speed - tau_E_k;
    }

    return 0;
}

/// Compute the squared speed after applying ground friction.
///
/// This function behaves similar to fric_speed(), except it operates on the
/// squared speed, making it more convenient to use if you are working on squared
/// speeds for performance reasons (e.g. by relying on fme_maxaccel_speed_C()). Compared
/// to fric_speed(), this function appears to be less performant when viewed in isolation.
/// However, if you are working on squared speeds, this function avoids computing
/// the square root for the most common case of geometric friction. Computing the
/// arithmetic friction requires a square root.
double fric_speedsq(double speedsq, double E, double tau_k)
{
    if (speedsq >= E * E) {
        const double tmp = 1 - tau_k;
        return speedsq * tmp * tmp;
    }

    const double tau_E_k = tau_k * E;
    const double tau_E_k_sq = tau_E_k * tau_E_k;
    if (speedsq >= tau_E_k_sq && speedsq >= 0.01) {
        return speedsq - 2 * std::sqrt(speedsq) * tau_E_k + tau_E_k_sq;
    }

    return 0;
}

/// Compute the speed after applying the FME.
///
/// This function runs at constant time.
double fme_speed(double speed, double costheta, double L, double ke_tau_M_A)
{
    double gamma2 = L - speed * costheta;
    if (gamma2 <= 0) {
        return speed;
    }

    double mu = ke_tau_M_A;
    if (gamma2 < mu) {
        mu = gamma2;
    }

    return std::sqrt(speed * (speed + 2 * mu * costheta) + mu * mu);
}

/// Compute the speed after applying the FME at maximum acceleration.
///
/// This function is easier to use than fme_maxaccel_speed_C(), but much less
/// performant for the common cases of zeta and 90 degrees strafing, due to
/// the need to compute a square root.
double fme_maxaccel_speed(double speed, double L, double ke_tau_M_A)
{
    if (ke_tau_M_A >= 0) {
        if (L <= ke_tau_M_A) {
            if (L >= 0) {
                return std::sqrt(speed * speed + L * L);
            }
            return speed;
        }

        const double tmp = L - ke_tau_M_A;
        if (tmp <= speed) {
            return std::sqrt(speed * speed + ke_tau_M_A * (L + tmp));
        }
        return speed + ke_tau_M_A;
    }

    if (-L < speed) {
        return speed - ke_tau_M_A;
    }
    return speed;
}

/// Compute a constant for adding to the squared speed at maximum acceleration.
///
/// The new squared speed after applying the FME can always be written in the form
/// of
///
/// \f[
///   \lVert\mathbf{v}'\rVert^2 = \lVert\mathbf{v}\rVert^2 + C(\lVert\mathbf{v}\rVert^2)
/// \f]
///
/// where \f$C\f$ is a function of the current squared speed in general. Though in
/// the zeta and the 90 degrees cases, \f$C\f$ is independent of the current squared
/// speed.
///
/// This function provides a significant performance boost over fme_maxaccel_speed() by
/// avoiding the need to compute the relatively expensive square root for the most common
/// cases of zeta and 90 degrees, at the expense of needing a square root when computing
/// the much less common linear and backward-linear cases. Since \f$C\f$ is independent
/// of the current speed in the zeta and 90 degrees cases, the value only needs to be
/// computed once in those cases, assuming the caller is certain that the
/// strafing type will never change. Subsequently, each iteration of the FME amounts only to
/// a single addition of \f$C\f$ to the current squared speed.
///
/// It is the responsibility of the user to use this function correctly. If you are unsure
/// how to use this function, opt for the more straightforward fme_maxaccel_speed() instead.
double fme_maxaccel_speed_C(double speedsq, double L, double ke_tau_M_A)
{
    if (ke_tau_M_A >= 0) {
        if (L <= ke_tau_M_A) {
            if (L >= 0) {
                return L * L;
            }
            return 0;
        }

        const double tmp = L - ke_tau_M_A;
        if (tmp * tmp <= speedsq) {
            return ke_tau_M_A * (L + tmp);
        }
        return (2 * std::sqrt(speedsq) + ke_tau_M_A) * ke_tau_M_A;
    }

    if (L >= 0 || L * L < speedsq) {
        return (ke_tau_M_A - 2 * std::sqrt(speedsq)) * ke_tau_M_A;
    }
    return 0;
}

// /// Compute the speed after applying the FME at minimum acceleration.
// ///
// double fme_minaccel_speed(double speed, double L, double ke_tau_M_A)
// {
//     if (ke_tau_M_A >= 0) {
//         if (L >= 0) {
//             if (L >= ke_tau_M_A) {
//                 return speed - ke_tau_M_A;
//             }

//             const double tmp = L - ke_tau_M_A;
//             if (L <= speed || -tmp <= speed) {
//                 return speed - ke_tau_M_A;
//             }
//             return L;
//         }

//         if (-L < speed) {
//             // TODO:
//         }
//     }


// }
