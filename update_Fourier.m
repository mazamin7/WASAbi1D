function p_next = update_Fourier(p_curr, p_prev, c, dt, dh, force, isDamped, alpha_abs)
% Computes p_next given p_curr, p_prev, c, dt, and dh using the formula
% p_next = 2 * p_curr - p_prev + (c * dt / dh)^2 * A * p_curr
%
% Inputs:
%   - p_curr: the current pressure values (a column vector)
%   - p_prev: the previous pressure values (a column vector)
%   - c: the speed of sound in the medium (a scalar)
%   - dt: the time step (a scalar)
%   - dh: the spatial step (a scalar)
%   - force: the applied force (a column vector)
%
% Output:
%   - p_next: the next pressure values (a column vector)

    N = length(p_curr);

    % modal analysis
    w = N/2;
    lx2 = w * w * dh * dh;
    
    cwt = zeros(w, 1);
    w2 = zeros(w, 1);
    
    for i = 1 : w
        ww = c * pi * sqrt(i^2 / lx2);
        w2(i) = ww^2;
        cwt(i) = cos(ww * dt);
    end

    p_prev_dct = dct(p_prev);
    p_curr_dct = dct(p_curr);
    p_next_dct = zeros(1,N);
    force_dct = dct(force);

    for i = 1 : w
        if isDamped == false
            p_next_dct(i) = 2.0 * p_curr_dct(i) * cwt(i) - p_prev_dct(i) ...
                + (2.0 * force_dct(i) / w2(i) ) * (1.0 - cwt(i));
        else
            lambda = w2(i);
            p_next_dct(i) = (2 - lambda * dt^2)/(1 + alpha_abs*dt/2) ...
                * p_curr_dct(i) - (1 - alpha_abs*dt/2) ...
                / (1 + alpha_abs*dt/2) * p_prev_dct(i);
        end
    end

    p_next = idct(p_next_dct);

end