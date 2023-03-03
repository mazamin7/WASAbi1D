function [p_next, q_next_dct] = update_Fourier(Fourier_data, p_curr, p_prev, force, q_curr_dct)
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
    isDamped = Fourier_data.isDamped;
    swt = Fourier_data.swt;
    alpha_abs = Fourier_data.alpha_abs;
    alpha2 = Fourier_data.alpha2;
    eatm = Fourier_data.eatm;
    w = Fourier_data.w;
    inv_w = Fourier_data.inv_w;
    inv_w2 = Fourier_data.inv_w2;

    DCT_type = 2;

    p_prev_dct = dct(p_prev,'Type',DCT_type);
    p_curr_dct = dct(p_curr,'Type',DCT_type);
    p_next_dct = zeros(1,N);
    force_dct = dct(force);
    % q_curr_dct = zeros(1,N);
    q_next_dct = zeros(1,N);

    for n = 1 : N
        if isDamped == false
            p_next_dct(n) = 2.0 * p_curr_dct(n) * cwt(n) - p_prev_dct(n) ...
                + (2.0 * force_dct(n) / w2(n) ) * (1.0 - cwt(n));
        else
            xe = force_dct(n) * inv_w2(n);
            p_next_dct(n) = xe + eatm(n) * ((p_curr_dct(n) - xe) * (cwt(n) + alpha_abs(n) * inv_w(n) * swt(n)) + swt(n) * inv_w(n) * q_curr_dct(n));
            q_next_dct(n) = eatm(n) * (q_curr_dct(n) * (cwt(n) - alpha_abs(n) * inv_w(n) * swt(n)) - (w(n) + alpha2(n) * inv_w(n)) * (p_curr_dct(n) - xe) * swt(n));
        end
    end

    p_next = idct(p_next_dct,'Type',DCT_type);

end