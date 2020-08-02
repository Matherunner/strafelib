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
