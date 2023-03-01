function p_next = update_FDTD(FDTD_data, p_curr, p_prev, force)
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
    boundCond1 = FDTD_data.boundCond1;
    boundCond2 = FDTD_data.boundCond2;
    c = FDTD_data.c;
    dt = FDTD_data.dt;
    dh = FDTD_data.dh;
    alpha_abs = FDTD_data.alpha_abs;
    isDamped = FDTD_data.isDamped;
    isPML = FDTD_data.isPML;
    sigma = FDTD_data.sigma;

    % Compute p_next using the formula
    if isDamped == false && isPML == false
        p_next = 2 * p_curr - p_prev + (c * dt / dh)^2 * A * p_curr + dt^2 * force;
    elseif isDamped == true
        p_next = (2 * p_curr - p_prev + alpha_abs*dt/2 * p_prev ...
            + (c * dt / dh)^2 * A * p_curr + dt^2 * force)/(1 + alpha_abs*dt/2);
    else
        p_next = 1 ./ (1 + dt * sigma) .* (2 * p_curr - p_prev + (c * dt / dh)^2 * A * p_curr ...
            + dt * sigma .* p_prev - dt * dt * sigma .* sigma .* p_curr);
    end

    if strcmp(boundCond1, "D")
        p_next(1) = 0;
    end

    if strcmp(boundCond2, "D")
        p_next(N) = 0;
    end

end