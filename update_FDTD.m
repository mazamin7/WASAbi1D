function p_next = update_FDTD(p_curr, p_prev, c, dt, dh, force, isDamped, alpha_abs, isBorrel)
% Computes p_next given p_curr, p_prev, c, dt, and dh using the formula
% p_next = 2 * p_curr - p_prev + (c * dt / dh)^2 * A * p_curr
%
% Inputs:
%   - p_curr: the current pressure values (a column vector)
%   - p_prev: the previous pressure values (a column vector)
%   - c: the speed of sound in the medium (a scalar)
%   - dt: the time step (a scalar)
%   - dh: the spatial step (a scalar)
%   - force: the applied force (a column vector)
%
% Output:
%   - p_next: the next pressure values (a column vector)

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;

    N = length(p_curr);
    
    A = sparse(N,N);
    
    if isBorrel == false
        A(1,1:4) = [delta + gamma, gamma + beta, beta + alpha, alpha];
        A(2,1:5) = [gamma + beta, delta + alpha, gamma, beta, alpha];
        A(3,1:6) = [beta + alpha, gamma, delta, gamma, beta, alpha];
        A(N-2,N-5:N) = [alpha, beta, gamma, delta, gamma, beta + alpha];
        A(N-1,N-4:N) = [alpha, beta, gamma, delta + alpha, gamma + beta];
        A(N,N-3:N) = [alpha, beta + alpha, gamma + beta, delta + gamma];
    else
        A(1,1:4) = [delta, 2*gamma, 2*beta, 2*alpha];
        A(2,1:5) = [gamma, delta + beta, gamma + alpha, beta, alpha];
        A(3,1:6) = [beta, gamma + alpha, delta, gamma, beta, alpha];
        A(N-2,N-5:N) = [alpha, beta, gamma, delta, gamma + alpha, beta];
        A(N-1,N-4:N) = [alpha, beta, gamma + alpha, delta + beta, gamma];
        A(N,N-3:N) = [2*alpha, 2*beta, 2*gamma, delta];
    end

    for i = 4:N-3
        A(i,i-3:i+3) = [alpha, beta, gamma, delta, gamma, beta, alpha];
    end

    % Compute p_next using the formula
    if isDamped == false
        p_next = 2 * p_curr - p_prev + (c * dt / dh)^2 * A * p_curr + dt^2 * force;
    else
        p_next = (2 * p_curr - p_prev + alpha_abs*dt/2 * p_prev ...
            + (c * dt / dh)^2 * A * p_curr + dt^2 * force)/(1 + alpha_abs*dt/2);
    end

end