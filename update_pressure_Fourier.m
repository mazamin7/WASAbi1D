function p_next = update_pressure_Fourier(Fourier_data, p_curr, p_prev, force, v_curr, override)
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
    order = Fourier_data.order;
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

    v_curr_dct = dct(v_curr,'Type',DCT_type);

	n = 2:N;

    % update solution in Fourier domain
    if order == 2 && override == false
        p_next_dct(n) = 2 * p_curr_dct(n) .* cwt(n) - p_prev_dct(n) ...
            + (2 * force_dct(n) ./ w2(n) ) .* (1 - cwt(n));
    else
        % the simulation code handles the correct instant of the force
        xe = force_dct(n) .* inv_w2(n);
        p_next_dct(n) = xe + eatm * ((p_curr_dct(n) - xe) .* (cwt(n) + alpha_abs * inv_w(n) .* swt(n)) + swt(n) .* inv_w(n) .* v_curr_dct(n));
    end

    n = 1;

    if order == 2 && override == false
        p_next_dct(n) = 2 * p_curr_dct(n) - p_prev_dct(n) + dt*dt * force_dct(n);
    else
        % the simulation code handles the correct instant of the force
        p_next_dct(n) = p_curr_dct(n) + dt * v_curr_dct(n);
    end

    % perform IDCT
    p_next = idct(p_next_dct,'Type',DCT_type);

end