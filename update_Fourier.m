function [p_next, p_prev_dct] = update_Fourier(Fourier_data, p_curr, p_prev, force)
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

    N = Fourier_data.N;
    w2 = Fourier_data.w2;
    cwt = Fourier_data.cwt;
    dt = Fourier_data.dt;
    alpha_abs = Fourier_data.alpha_abs;
    isDamped = Fourier_data.isDamped;

    p_prev_dct = dct(p_prev);
    p_curr_dct = dct(p_curr);
    p_next_dct = zeros(1,N);
    force_dct = dct(force);

    for n = 1 : N
        if isDamped == false
            p_next_dct(n) = 2.0 * p_curr_dct(n) * cwt(n) - p_prev_dct(n) ...
                + (2.0 * force_dct(n) / w2(n) ) * (1.0 - cwt(n));
        else
            lambda = w2(n);
            p_next_dct(n) = (2 - lambda * dt^2)/(1 + alpha_abs*dt/2) ...
                * p_curr_dct(n) - (1 - alpha_abs*dt/2) ...
                / (1 + alpha_abs*dt/2) * p_prev_dct(n);
        end
    end

    p_next = idct(p_next_dct);

end