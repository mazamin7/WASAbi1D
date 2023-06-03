function [stable] = check_stability(c, dt, dh, order, fourier, diss, DD)

    % check stability
    if fourier == false % FDTD
        stable = check_stability_FDTD(c, dt, dh, order, diss, true);
    elseif fourier == true && DD == false % Fourier no DD
        stable = check_stability_Fourier(diss);
    elseif fourier == true && DD == true % Fourier DD
        cond1 = check_stability_Fourier(diss);
        cond2 = check_stability_FDTD(c, dt, dh, order, diss, false);
        stable = cond1 && cond2;
    end

end