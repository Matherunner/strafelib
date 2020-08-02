/** @file */

#include <cmath>

constexpr const double JUMP_SPEED = 268.3281572999748;

/// Compute the speed after applying ground friction.
///
/// This function runs at constant time. The caller is responsible of computing
/// \p tau_k correctly by considering the values of the entity friction and
/// the edgefriction.
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
