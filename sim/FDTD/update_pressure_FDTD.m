function p_next = update_pressure_FDTD(data, p_curr, p_prev, force, v_curr, override)
% Computes p_next given p_curr, p_prev, force and FDTD_data
%
% Inputs:
%   - FDTD_data: data structure containing data needed for the
%   computation
%   - p_curr: the current pressure values (a column vector)
%   - p_prev: the previous pressure values (a column vector)
%   - force: the applied force (a column vector)
%
% Output:
%   - p_next: the next pressure values (a column vector)

    N = data.N;
    A = data.A;
    c = data.c;
    dt = data.dt;
    dh = data.dh;
    alpha_abs = data.alpha_abs;
    temp_order = data.temp_order;

    if temp_order == 2 && override == false
        % Compute p_next using the formula
        % current force
        p_next = (2 * p_curr - (1 - dt*alpha_abs) * p_prev ...
            + (c * dt / dh)^2 * A * p_curr + dt^2 * force)/(1 + dt*alpha_abs);
    else
        % no force
        p_next = p_curr + dt * v_curr;
    end

end