function [p_next, v_next] = update_FDTD(data, p_curr, p_prev, force, v_curr, v_prev, g1, g2)
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
    bc_left = data.bc_left;
    bc_right = data.bc_right;
    c = data.c;
    dt = data.dt;
    dh = data.dh;
    alpha_abs = data.alpha_abs;
    isPML = data.isPML;
    sigma = data.sigma;
    order = data.order;

    % Extending solutions to include ghost points
    p_curr_old = p_curr;
    p_curr = zeros(N+4,1);
    p_curr(3:end-2) = p_curr_old;

    p_prev_old = p_prev;
    p_prev = zeros(N+4,1);
    p_prev(3:end-2) = p_prev_old;

    q_curr_old = v_curr;
    v_curr = zeros(N+4,1);
    v_curr(3:end-2) = q_curr_old;

    q_prev_old = v_prev;
    v_prev = zeros(N+4,1);
    v_prev(3:end-2) = q_prev_old;

    force_old = force';
    force = zeros(N+4,1);
    force(3:end-2) = force_old;

    sigma_old = sigma';
    sigma = zeros(N+4,1);
    sigma(3:end-2) = sigma_old;

    % Imposing b.c. on the left
    if strcmp(bc_left, "D")
        p_curr(1:3) = [1 1 1] * g1;
    elseif strcmp(bc_left, "N")
        p_curr(3) = p_curr(3) - 0.5 * dh * g1;
        p_curr(4) = p_curr(4) - 1.5 * dh * g1;
        p_curr(5) = p_curr(5) - 2.5 * dh * g1;
    end

    % Imposing b.c. on the right
    if strcmp(bc_right, "D")
        p_curr(end-2:end) = [1 1 1] * g2;
    elseif strcmp(bc_right, "N")
        p_curr(end-5) = p_curr(end-5) - 0.5 * dh * g2;
        p_curr(end-4) = p_curr(end-4) - 1.5 * dh * g2;
        p_curr(end-3) = p_curr(end-3) - 2.5 * dh * g2;
    end

    if order == 2
        % Compute p_next using the formula
        if isPML == false
            p_next = (2 * p_curr - p_prev + alpha_abs*dt * p_prev ...
                + (c * dt / dh)^2 * A * p_curr + dt^2 * force)/(1 + alpha_abs*dt);
        else % isPML == true
            p_next = (2 * p_curr - p_prev + (c * dt / dh)^2 * A * p_curr ...
                + dt * sigma .* p_prev - dt^2 * sigma^2 .* p_curr) ./ (1 + dt * sigma);
        end
    
        % we will update after post-merge
        % we could use: v_next = (p_next - p_curr)/dt;
        v_next = v_curr * 0;
    elseif order == 1
        % Compute p_next and q_next using the formula
%         p_next = 2 * dt * v_curr + p_prev;
%         v_next = 2 * c^2 * dt / dh^2 * A * p_curr + v_prev - 4 * dt * alpha_abs * v_curr + 2 * dt * force;

        % new proposal (bad)
%         v_next = (1-2*alpha_abs*dt) * v_curr + c^2 * dt / dh^2 * A * p_curr;
%         p_next = dt * v_curr + p_curr;

        % v^{n+1} = (v^{n-1} + 2dt*c^2/dh^2*p^{n} + 2dt*F^{n+1})/(1+4dt*alpha)
        % I SHOULD USE NEXT FORCE
        p_next = 2 * dt * v_curr + p_prev;
        v_next = (v_prev + 2 * c^2 * dt / dh^2 * A * p_curr + 2 * dt * force)/(1 + 4*dt*alpha_abs);
    end

    % Truncating ghost points
    p_next = p_next(3:end-2);
    v_next = v_next(3:end-2);

end