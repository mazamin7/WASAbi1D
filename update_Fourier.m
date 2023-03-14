function [p_next, q_next_dct] = update_Fourier(Fourier_data, p_curr, p_prev, force, q_curr_dct, g)
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
    boundCond1 = Fourier_data.boundCond1;
    boundCond2 = Fourier_data.boundCond2;
    isLeft = Fourier_data.isLeft;
    dh = Fourier_data.dh;
    c = Fourier_data.c;
    C1 = Fourier_data.C1;
    C2 = Fourier_data.C2;

    % performing DCT
    DCT_type = 2;

    p_prev_dct = dct(p_prev,'Type',DCT_type);
    p_curr_dct = dct(p_curr,'Type',DCT_type);
    p_next_dct = zeros(N,1);
    force_dct = dct(force);
    % q_curr_dct = zeros(N,1);
    q_next_dct = zeros(N,1);

    % update solution in Fourier domain
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

    % perform IDCT
    p_next = idct(p_next_dct,'Type',DCT_type);

    % impose boundary conditions
    if strcmp(boundCond1, "D")
        % reversing Neumann
        p_next = p_next + (c * dt / dh)^2 * C1 * p_curr;

        % imposing Dirichlet
        p_next(1) = g;
    elseif isLeft && strcmp(boundCond1, "N")
        % imposing Neumann deviation from homogeneity
        p_next(1) = p_next(1) - dh * 0.5 * g;
        p_next(2) = p_next(2) - dh * 1.5 * g;
        p_next(3) = p_next(3) - dh * 2.5 * g;
    end

    if strcmp(boundCond2, "D")
        % reversing Neumann
        p_next = p_next + (c * dt / dh)^2 * C2 * p_curr;

        % imposing Dirichlet
        p_next(N) = g;
    elseif isLeft == false && strcmp(boundCond2, "N")
        % imposing Neumann deviation from homogeneity
        p_next(N-2) = p_next(N-2) + dh * 2.5 * g;
        p_next(N-1) = p_next(N-1) + dh * 1.5 * g;
        p_next(N) = p_next(N) + dh * 0.5 * g;
    end

end