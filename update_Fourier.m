function [p_next, q_next_dct] = update_Fourier(Fourier_data, p_curr, p_prev, force, q_curr_dct, exact_damped)
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

    % retrieving data
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

    % performing DCT
    DCT_type = 2;

    p_prev_dct = dct(p_prev,'Type',DCT_type);
    p_curr_dct = dct(p_curr,'Type',DCT_type);
    p_next_dct = zeros(N,1);
    force_dct = dct(force);
    % q_curr_dct = zeros(N,1);
    q_next_dct = zeros(N,1);

    % update solution in Fourier domain
    for n = 2 : N
        if isDamped == false
            p_next_dct(n) = 2 * p_curr_dct(n) * cwt(n) - p_prev_dct(n) ...
                + (2 * force_dct(n) / w2(n) ) * (1 - cwt(n));
        elseif isDamped == true && exact_damped == false
            p_next_dct(n) = (2 - w2(n) * dt^2)/(1 + alpha_abs*dt/2) ...
                * p_curr_dct(n) - (1 - alpha_abs*dt/2) ...
                / (1 + alpha_abs*dt/2) * p_prev_dct(n);
        elseif isDamped == true && exact_damped == true
            xe = force_dct(n) * inv_w2(n);
            p_next_dct(n) = xe + eatm * ((p_curr_dct(n) - xe) * (cwt(n) + alpha_abs * inv_w(n) * swt(n)) + swt(n) * inv_w(n) * q_curr_dct(n));
            q_next_dct(n) = eatm * (q_curr_dct(n) * (cwt(n) - alpha_abs * inv_w(n) * swt(n)) - (w(n) + alpha2 * inv_w(n)) * (p_curr_dct(n) - xe) * swt(n));
        end
    end

    assert(p_next_dct(1) == 0);

    % perform IDCT
    p_next = idct(p_next_dct,'Type',DCT_type);

end