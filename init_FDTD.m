function FDTD_data = init_FDTD(N, c, dt, dh, isDamped, alpha_abs, isBorrel, boundCond1, boundCond2, isPML)

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;
    
    A = sparse(N+4,N+4);

    for i = 4:N+1
        A(i,i-3:i+3) = [alpha, beta, gamma, delta, gamma, beta, alpha];
    end

    if strcmp(boundCond1, "N")
        if isBorrel == false || (isBorrel == true && isLeft)
            A(3,3:6) = [delta + gamma, gamma + beta, beta + alpha, alpha];
            A(4,3:7) = [gamma + beta, delta + alpha, gamma, beta, alpha];
            A(5,3:8) = [beta + alpha, gamma, delta, gamma, beta, alpha];
        else
            A(3,3:6) = [delta, 2*gamma, 2*beta, 2*alpha];
            A(4,3:7) = [gamma, delta + beta, gamma + alpha, beta, alpha];
            A(5,3:8) = [beta, gamma + alpha, delta, gamma, beta, alpha];
        end
    end

    if strcmp(boundCond2, "N")
        if isBorrel == false || (isBorrel == true && isLeft == false)
            A(end-4,N-3:N+2) = [alpha, beta, gamma, delta, gamma, beta + alpha];
            A(end-3,N-2:N+2) = [alpha, beta, gamma, delta + alpha, gamma + beta];
            A(end-2,N-1:N+2) = [alpha, beta + alpha, gamma + beta, delta + gamma];
        else
            A(end-4,N-3:N+2) = [alpha, beta, gamma, delta, gamma + alpha, beta];
            A(end-3,N-2:N+2) = [alpha, beta, gamma + alpha, delta + beta, gamma];
            A(end-2,N-1:N+2) = [2*alpha, 2*beta, 2*gamma, delta];
        end
    end

    sigma = zeros(N,1);
    
    for i = 2:N
        sigma(i) = sigma(i-1) + 2;
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
    FDTD_data.isBorrel = isBorrel;
    FDTD_data.isPML = isPML;
    FDTD_data.sigma = sigma;

end