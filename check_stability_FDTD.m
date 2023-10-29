function [stable] = check_stability_FDTD(c, dt, dh, space_order, diss, asymptotic)
%
% Only works for space_order == 6
%

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;

    var_1511 = 1/2 * (abs(alpha) + abs(beta) + abs(gamma) - delta/2);

    CFL = 1 / sqrt(var_1511);

    if asymptotic
        stable = dt < dh / c * CFL; % CFL condition
    else
        stable = dt <= dh / c * CFL;
    end

    stable = stable && diss < 1;

end