function [stable] = check_stability(c, dt, dh, alpha_abs, order, xi, nu, fourier, nu_fourier, DD)

    % check stability
    if fourier == false % FDTD
        stable = check_stability_FDTD(c, dt, dh, alpha_abs, order, xi, nu, true);
    elseif fourier == true && DD == false % Fourier no DD
        stable = check_stability_Fourier(nu_fourier);
    elseif fourier == true && DD == true % Fourier DD
        cond1 = check_stability_Fourier(nu_fourier);
        cond2 = check_stability_FDTD(c, dt, dh, alpha_abs, order, 1, nu_fourier, false);
        stable = cond1 && cond2;
    end

end