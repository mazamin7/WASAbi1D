function Fourier_data = init_Fourier(N, c, dt, dh, isDamped, alpha_abs)

    % modal analysis
    w = N/2;
    lx2 = w * w * dh * dh;
    
    cwt = zeros(N, 1);
    w2 = zeros(N, 1);
    
    for i = 1 : w
        ww = c * pi * sqrt(i^2 / lx2);
        w2(i) = ww^2;
        cwt(i) = cos(ww * dt);
    end

    Fourier_data.N = N;
    Fourier_data.w2 = w2;
    Fourier_data.cwt = cwt;
    Fourier_data.dt = dt;
    Fourier_data.alpha_abs = alpha_abs;
    Fourier_data.isDamped = isDamped;

end