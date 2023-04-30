function data = init_FDTD_1ord(len_x, c, dt, dh, alpha_abs)

    N = floor(len_x/dh);
    N = 2 * floor(N/2);

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;
    
    A = sparse(N,N);

    for i = 4:N-3
        A(i,i-3:i+3) = [alpha, beta, gamma, delta, gamma, beta, alpha];
    end

    A(1,1:4) = [delta + gamma, gamma + beta, beta + alpha, alpha];
    A(2,1:5) = [gamma + beta, delta + alpha, gamma, beta, alpha];
    A(3,1:6) = [beta + alpha, gamma, delta, gamma, beta, alpha];
    A(end-2,N-5:N) = [alpha, beta, gamma, delta, gamma, beta + alpha];
    A(end-1,N-4:N) = [alpha, beta, gamma, delta + alpha, gamma + beta];
    A(end,N-3:N) = [alpha, beta + alpha, gamma + beta, delta + gamma];

    data.N = N;
    data.A = A;
    data.c = c;
    data.dt = dt;
    data.dh = dh;
    data.alpha_abs = alpha_abs;

end