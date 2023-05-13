function [stable] = check_stability(len_x, c, dt, dh, alpha_abs, order, diss)

    N = floor(len_x/dh);

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;
    
    K = sparse(N,N);

    for i = 4:N-3
        K(i,i-3:i+3) = [alpha, beta, gamma, delta, gamma, beta, alpha];
    end

    K(1,1:4) = [delta + gamma, gamma + beta, beta + alpha, alpha];
    K(2,1:5) = [gamma + beta, delta + alpha, gamma, beta, alpha];
    K(3,1:6) = [beta + alpha, gamma, delta, gamma, beta, alpha];

    K(end-2,N-5:N) = [alpha, beta, gamma, delta, gamma, beta + alpha];
    K(end-1,N-4:N) = [alpha, beta, gamma, delta + alpha, gamma + beta];
    K(end,N-3:N) = [alpha, beta + alpha, gamma + beta, delta + gamma];

    % check stability
    if order == 1
        id_mat = diag(ones(N,1));
        zero_mat = zeros(N,N);

        K_overline = [-2*alpha_abs*id_mat, c^2/dh^2*K;
                      id_mat, zero_mat];

        id_mat = eye(2*N);
        zero_mat = zeros(2*N,2*N);
    
        B = [2*dt*K_overline, id_mat * diss;
             id_mat, zero_mat];
    
        BD = eigs(B, 4*N);
        BD = max(abs(BD));
    
        stable = BD < 1; % check that largest abs is less than 1
    else
        var_1511 = 1/2 * (abs(alpha) + abs(beta) + abs(gamma) - delta/2);
        CFL = 1 / sqrt(var_1511);
        stable = dt < dh / c * CFL; % CFL condition
    end

end