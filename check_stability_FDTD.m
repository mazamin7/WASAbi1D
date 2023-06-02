function [stable] = check_stability_FDTD(c, dt, dh, alpha_abs, order, xi, nu, asymptotic)

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;

    var_1511 = 1/2 * (abs(alpha) + abs(beta) + abs(gamma) - delta/2);

    if order == 1
        CFL = 0.5 / sqrt(var_1511);
    elseif order == 2
        CFL = 1 / sqrt(var_1511);
    end

    if asymptotic
        stable = dt < dh / c * CFL; % CFL condition
    else
        stable = dt <= dh / c * CFL;
    end

end