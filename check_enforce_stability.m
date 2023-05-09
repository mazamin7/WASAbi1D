function [stable, redux] = check_enforce_stability(len_x, c, dt, dh, alpha_abs, order, diss)

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

        K_overline = [-4*dt*alpha_abs*id_mat, 2*dt*c^2/dh^2*K;
            2*dt*id_mat, zero_mat];

        max_eig_K_overline = max(abs(eigs(K_overline,2*N)));

        if max_eig_K_overline < 1
            id_mat = diag(ones(2*N,1));
            zero_mat = zeros(2*N,2*N);
            B = [K_overline, id_mat;
                id_mat, zero_mat];
    
            BD = eigs(B, 4*N);
            disp(max(abs(BD)))
            % forcing the scheme to converge
            redux = diss / max(abs(BD));
    
            B = [K_overline, id_mat * redux;
                id_mat, zero_mat];
    
            BD = eigs(B, 4*N);
            disp(abs(BD)) % get abs of eigenvalues
            disp(max(abs(BD))) % get max abs of eigenvalues
    
            stable = max(abs(BD)) < 1; % check that largest abs is less than 1
        else
            % dt too small, can't converge
            stable = 0;
            redux = 0;
        end
    else
        stable = dt < dh / c; % CFL condition
        redux = 0;
    end

end