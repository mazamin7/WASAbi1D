function FDTD_data = init_FDTD(N, c, dt, dh, isDamped, alpha_abs, isBorrel, boundCond1, boundCond2)

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;
    
    A = sparse(N,N);
    
    if strcmp(boundCond1, "N")
        if isBorrel == false
            A(1,1:4) = [delta + gamma, gamma + beta, beta + alpha, alpha];
            A(2,1:5) = [gamma + beta, delta + alpha, gamma, beta, alpha];
            A(3,1:6) = [beta + alpha, gamma, delta, gamma, beta, alpha];
        else
            A(1,1:4) = [delta, 2*gamma, 2*beta, 2*alpha];
            A(2,1:5) = [gamma, delta + beta, gamma + alpha, beta, alpha];
            A(3,1:6) = [beta, gamma + alpha, delta, gamma, beta, alpha];
        end
    elseif strcmp(boundCond1, "D")
        A(1,1:4) = [delta, gamma, beta, alpha];
        A(2,1:5) = [gamma, delta, gamma, beta, alpha];
        A(3,1:6) = [beta, gamma, delta, gamma, beta, alpha];
    end

    if strcmp(boundCond2, "N")
        if isBorrel == false
            A(N-2,N-5:N) = [alpha, beta, gamma, delta, gamma, beta + alpha];
            A(N-1,N-4:N) = [alpha, beta, gamma, delta + alpha, gamma + beta];
            A(N,N-3:N) = [alpha, beta + alpha, gamma + beta, delta + gamma];
        else
            A(N-2,N-5:N) = [alpha, beta, gamma, delta, gamma + alpha, beta];
            A(N-1,N-4:N) = [alpha, beta, gamma + alpha, delta + beta, gamma];
            A(N,N-3:N) = [2*alpha, 2*beta, 2*gamma, delta];
        end
    elseif strcmp(boundCond2, "D")
        A(N-2,N-5:N) = [alpha, beta, gamma, delta, gamma, beta];
        A(N-1,N-4:N) = [alpha, beta, gamma, delta, gamma];
        A(N,N-3:N) = [alpha, beta, gamma, delta];
    end

    for i = 4:N-3
        A(i,i-3:i+3) = [alpha, beta, gamma, delta, gamma, beta, alpha];
    end

    FDTD_data.N = N;
    FDTD_data.A = A;
    FDTD_data.boundCond1 = boundCond1;
    FDTD_data.boundCond2 = boundCond2;
    FDTD_data.c = c;
    FDTD_data.dt = dt;
    FDTD_data.dh = dh;
    FDTD_data.alpha_abs = alpha_abs;
    FDTD_data.isDamped = isDamped;

end