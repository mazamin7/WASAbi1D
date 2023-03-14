function Fourier_data = init_Fourier(N, c, dt, dh, isDamped, alpha, boundCond, isLeft)

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;

    % modal analysis
    inv_lx = 1 / ((N - 1) * dh);
    
    cwt = zeros(N, 1);
    w2 = zeros(N, 1);
    swt = zeros(N, 1);
    alpha_abs = zeros(N, 1);
    alpha2 = zeros(N, 1);
    eatm = zeros(N, 1);
    w = zeros(N, 1);
    inv_w = zeros(N, 1);
    inv_w2 = zeros(N, 1);
    
    for n = 1 : N
        ww = n * pi * c * inv_lx;
        % disp(['mode ' num2str(n) ', freq [rad/s] = ' num2str(ww)])
        w(n) = ww;
        inv_w(n) = 1/ww;
        inv_w2(n) = inv_w(n) * inv_w(n);
        w2(n) = ww * ww;
        cwt(n) = cos(ww * dt);
        swt(n) = sin(ww * dt);
        alpha_abs(n) = alpha;
        alpha2(n) = alpha_abs(n) * alpha_abs(n);
        eatm(n) = exp(- alpha_abs(n) * dt);
    end

    if isLeft
        boundCond1 = boundCond;
        boundCond2 = "N";
    else
        boundCond1 = "N";
        boundCond2 = boundCond;
    end

    C1 = sparse(N,N);
    C1(1,1:3) = -[gamma beta alpha];
    C1(2,1:2) = -[beta alpha];
    C1(3,1) = -alpha;
    
    C2 = sparse(N,N);
    C2(N-2,N) = -alpha;
    C2(N-1,N-1:N) = -[beta alpha];
    C2(N,N-2:N) = -[gamma beta alpha];

    Fourier_data.N = N;
    Fourier_data.w2 = w2;
    Fourier_data.cwt = cwt;
    Fourier_data.dt = dt;
    Fourier_data.isDamped = isDamped;
    Fourier_data.swt = swt;
    Fourier_data.alpha_abs = alpha_abs;
    Fourier_data.alpha2 = alpha2;
    Fourier_data.eatm = eatm;
    Fourier_data.w = w;
    Fourier_data.inv_w = inv_w;
    Fourier_data.inv_w2 = inv_w2;
    Fourier_data.boundCond1 = boundCond1;
    Fourier_data.boundCond2 = boundCond2;
    Fourier_data.isLeft = isLeft;
    Fourier_data.dh = dh;
    Fourier_data.C1 = C1;
    Fourier_data.C2 = C2;
    Fourier_data.c = c;

end