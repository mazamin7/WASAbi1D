function Fourier_data = init_Fourier(N, c, dt, dh, isDamped, alpha)

    % modal analysis
    inv_lx = 1 / ((N - 1) * dh);
    
    cwt = zeros(N, 1);
    w2 = zeros(N, 1);
    swt = zeros(N, 1);
    alpha_abs = 0;
    alpha2 = 0;
    eatm = 0;
    w = zeros(N, 1);
    inv_w = zeros(N, 1);
    inv_w2 = zeros(N, 1);

    alpha_abs = alpha;
    alpha2 = alpha_abs * alpha_abs;
    eatm = exp(- alpha_abs * dt);
    
    for n = 1 : N
        ww = (n-1) * pi * c * inv_lx;
        % disp(['mode ' num2str(n) ', freq [rad/s] = ' num2str(ww)])
        w(n) = ww;
        inv_w(n) = 1/ww;
        inv_w2(n) = inv_w(n) * inv_w(n);
        w2(n) = ww * ww;
        cwt(n) = cos(ww * dt);
        swt(n) = sin(ww * dt);
    end

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

end