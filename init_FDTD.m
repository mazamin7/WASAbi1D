function data = init_FDTD(len_x, c, dt, dh, alpha_abs, bc_left, bc_right, isPML, order)

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

    data.N = N;
    data.A = A;
    data.bc_left = bc_left;
    data.bc_right = bc_right;
    data.c = c;
    data.dt = dt;
    data.dh = dh;
    data.alpha_abs = alpha_abs;
    data.isPML = isPML;
    data.sigma = sigma;
    data.order = order;

end