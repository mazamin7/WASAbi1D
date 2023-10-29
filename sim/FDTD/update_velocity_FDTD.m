function v_next = update_velocity_FDTD(data, p_next, p_curr, p_prev, force_next, v_curr, override)
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
        % no force
        v_next = (p_next - p_curr)/dt;
    else
        % next force
        v_next = (v_curr + c^2 * dt / dh^2 * A * p_next + dt * force_next)/(1 + 2*dt*alpha_abs);
    end

end