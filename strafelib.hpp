/** @file */

#include <cmath>

/// The initial jumping speed, before gravity is applied.
constexpr const double JUMP_SPEED = 268.3281572999748;

/// Compute the dot product of two vectors.
///
template<int N>
inline double dot_product(const double *__restrict a, const double *__restrict b)
{
    double res = 0;
    for (int i = 0; i < N; ++i) {
        res += a[i] * b[i];
    }
    return res;
}

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

/// Compute the velocity after applying ground friction.
///
/// \p vel must have at least two elements.
///
/// The caller is responsible of ensuring \p speed matches the 2D norm
/// of \p vel. By giving you the responsibility of computing the speed,
/// this function avoids computing any square roots.
void fric_vel(double *__restrict vel, double speed, double E, double tau_k)
{
    if (speed >= E) {
        const double tmp = 1 - tau_k;
        vel[0] *= tmp;
        vel[1] *= tmp;
        return;
    }

    const double tau_E_k = tau_k * E;
    if (speed >= tau_E_k && speed >= 0.1) {
        const double tmp = tau_E_k / speed;
        vel[0] -= vel[0] * tmp;
        vel[1] -= vel[1] * tmp;
        return;
    }

    vel[0] = 0;
    vel[1] = 0;
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
    const double gamma2 = L - speed * costheta;
    if (gamma2 <= 0) {
        return speed;
    }

    double mu = ke_tau_M_A;
    if (gamma2 < mu) {
        mu = gamma2;
    }

    return std::sqrt(speed * (speed + 2 * mu * costheta) + mu * mu);
}

/// Compute the velocity after applying the FME.
///
/// The caller is responsible of ensuring \p costheta and \p sintheta are consistent.
/// \p speed must also match the 2D magnitude of \p vel. \p vel must have at least two elements.
///
/// The apparent redundancy of providing both \p vel and \p speed is to improve the
/// efficiency. Under most situations, the user of this function needs to compute a
/// square root to obtain the speed (the norm of velocity). This speed is needed by
/// a function that computes \f$\cos\theta\f$ and this function. To avoid computing
/// the square root twice, we leave the responsibility to the caller.
void fme_vel_theta(double *__restrict vel, double speed, double costheta, double sintheta, double L, double ke_tau_M_A)
{
    const double gamma2 = L - speed * costheta;
    if (gamma2 <= 0) {
        return;
    }

    const double mu = std::min(ke_tau_M_A, gamma2);
    const double tmp = mu / speed;
    const double ax = vel[0] * costheta + vel[1] * sintheta;
    const double ay = vel[1] * costheta - vel[0] * sintheta;
    vel[0] += tmp * ax;
    vel[1] += tmp * ay;
}

/// Compute the \f$\cos\theta\f$ for maximum acceleration.
///
/// There is no point accepting a squared speed, because the very common
/// zeta strafing case requires computing a square root anyway.
void fme_maxaccel_cossin_theta(double speed, double L, double ke_tau_M_A, double *__restrict costheta, double *__restrict sintheta)
{
    if (ke_tau_M_A >= 0) {
        if (L <= ke_tau_M_A) {
            if (L >= 0) {
                *costheta = 0;
                *sintheta = 1;
                return;
            }
            *costheta = 1;
            *sintheta = 0;
            return;
        }

        const double tmp = L - ke_tau_M_A;
        if (tmp <= speed) {
            const double ct = tmp / speed;
            *costheta = ct;
            *sintheta = std::sqrt(1 - ct * ct);
            return;
        }
        *costheta = 1;
        *sintheta = 0;
        return;
    }

    if (-L < speed) {
        *costheta = -1;
        *sintheta = 0;
        return;
    }
    *costheta = 1;
    *sintheta = 0;
    return;
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

/// Compute the speed after applying the FME at minimum acceleration.
///
double fme_minaccel_speed(double speed, double L, double ke_tau_M_A)
{
    if (ke_tau_M_A >= 0) {
        if (L >= 0) {
            if (L >= ke_tau_M_A) {
                return std::fabs(speed - ke_tau_M_A);
            }

            const double tmp = L - ke_tau_M_A;
            if (L <= speed || -tmp <= speed) {
                return std::fabs(speed - ke_tau_M_A);
            }
            return L;
        }

        if (-L < speed) {
            if (ke_tau_M_A - L <= speed) {
                return std::fabs(speed - ke_tau_M_A);
            }
            return L; 
        }
        return speed;
    }

    return speed;
}

/// Compute the velocity after colliding with a hyperplane.
///
/// The caller is responsible of ensuring \p n is a unit vector.
template<int N>
void collision_vel(double *__restrict v, const double *__restrict n, double b)
{
    const double tmp = b * dot_product<N>(v, n);
    for (int i = 0; i < N; ++i) {
        v[i] -= tmp * n[i];
    }
}

/// Compute the hunting velocity of a snark.
///
/// \p dir must be a unit vector.
template<int N>
void snark_hunt_vel(double *__restrict v, const double *__restrict dir)
{
    const double speed = std::sqrt(dot_product<N>(v, v));
    double tmp = 1.2;
    if (speed > 95. / 3) {
        tmp = 50 / (speed + 10);
    }
    for (int i = 0; i < N; ++i) {
        v[i] = tmp * v[i] + 300 * dir[i];
    }
}

/// Compute the player frame time given the game frame time.
///
inline double tau_g_to_p(double tau_g)
{
    return 0.001 * std::floor(1000 * tau_g);
}
