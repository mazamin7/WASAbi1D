function FDTD_data = init_FDTD(len_x, c, dt, dh, alpha_abs, bc_left, bc_right, isPML)

    N = floor(len_x/dh);
    N = 2 * floor(N/2);

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;
    
    A = sparse(N+4,N+4);

    for i = 4:N+1
        A(i,i-3:i+3) = [alpha, beta, gamma, delta, gamma, beta, alpha];
    end

    if strcmp(bc_left, "N")
        A(3,3:6) = [delta + gamma, gamma + beta, beta + alpha, alpha];
        A(4,3:7) = [gamma + beta, delta + alpha, gamma, beta, alpha];
        A(5,3:8) = [beta + alpha, gamma, delta, gamma, beta, alpha];
    end

    if strcmp(bc_right, "N")
        A(end-4,N-3:N+2) = [alpha, beta, gamma, delta, gamma, beta + alpha];
        A(end-3,N-2:N+2) = [alpha, beta, gamma, delta + alpha, gamma + beta];
        A(end-2,N-1:N+2) = [alpha, beta + alpha, gamma + beta, delta + gamma];
    end

    sigma = zeros(N,1);
    
    for i = 2:N
        sigma(i) = 2*i;
    end

    FDTD_data.N = N;
    FDTD_data.A = A;
    FDTD_data.bc_left = bc_left;
    FDTD_data.bc_right = bc_right;
    FDTD_data.c = c;
    FDTD_data.dt = dt;
    FDTD_data.dh = dh;
    FDTD_data.alpha_abs = alpha_abs;
    FDTD_data.isPML = isPML;
    FDTD_data.sigma = sigma;

end