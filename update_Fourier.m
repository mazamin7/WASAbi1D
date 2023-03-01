function [p_next, p_prev_dct] = update_Fourier(Fourier_data, p_curr, p_prev, force)
% Computes p_next given p_curr, p_prev, force and Fourier_data
%
% Inputs:
%   - Fourier_data: data structure containing data needed for the
%   computation
%   - p_curr: the current pressure values (a column vector)
%   - p_prev: the previous pressure values (a column vector)
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

    DCT_type = 2;

    p_prev_dct = dct(p_prev,'Type',DCT_type);
    p_curr_dct = dct(p_curr,'Type',DCT_type);
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

    p_next = idct(p_next_dct,'Type',DCT_type);

end