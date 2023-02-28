function Fourier_data = init_Fourier(N, c, dt, dh, isDamped, alpha_abs)

    % modal analysis
    inv_lx = 1 / ((N - 1) * dh);
    
    cwt = zeros(N, 1);
    w2 = zeros(N, 1);
    
    for n = 1 : N
        ww = n * pi * c * inv_lx;
        % disp(['mode ' num2str(n) ', freq [rad/s] = ' num2str(ww)])
        w2(n) = ww * ww;
        cwt(n) = cos(ww * dt);
    end

    Fourier_data.N = N;
    Fourier_data.w2 = w2;
    Fourier_data.cwt = cwt;
    Fourier_data.dt = dt;
    Fourier_data.alpha_abs = alpha_abs;
    Fourier_data.isDamped = isDamped;

end