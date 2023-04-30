function [p_next, q_next] = update_FDTD_1ord(FDTD_data, p_curr, p_prev, force, q_curr, q_prev)
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

    N = FDTD_data.N;
    A = FDTD_data.A;
    c = FDTD_data.c;
    dt = FDTD_data.dt;
    dh = FDTD_data.dh;
    alpha_abs = FDTD_data.alpha_abs;

    % Compute p_next and q_next using the formula
    q_next = 2 * c^2 * dt / dh^2 * A * p_curr + q_prev - 4 * dt * alpha_abs * q_curr + 2 * dt * force;
    p_next = 2 * dt * q_curr + p_prev;

end