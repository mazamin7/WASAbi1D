function Fourier_data = init_Fourier(len_x, c, dt, dh, temp_order, alpha_abs)

    % modal analysis
    inv_lx = 1 / len_x;

    N = floor(len_x/dh);
    
    cwt = zeros(N, 1);
    w2 = zeros(N, 1);
    swt = zeros(N, 1);
    w = zeros(N, 1);
    inv_w = zeros(N, 1);
    inv_w2 = zeros(N, 1);

    alpha2 = alpha_abs * alpha_abs;
    eatm = exp(- alpha_abs * dt);
    
	% first element is empty but no one cares
	n = 2:N;
	
    w_0 = pi * (n-1) * c * inv_lx;
    w(n) = w_0;
    inv_w(n) = 1./w(n);
    inv_w2(n) = inv_w(n) .* inv_w(n);
    w2(n) = w(n) .* w(n);
    cwt(n) = cos(w(n) * dt);
    swt(n) = sin(w(n) * dt);

    Fourier_data.N = N;
    Fourier_data.w2 = w2;
    Fourier_data.cwt = cwt;
    Fourier_data.dt = dt;
    Fourier_data.temp_order = temp_order;
    Fourier_data.swt = swt;
    Fourier_data.alpha_abs = alpha_abs;
    Fourier_data.alpha2 = alpha2;
    Fourier_data.eatm = eatm;
    Fourier_data.w = w;
    Fourier_data.inv_w = inv_w;
    Fourier_data.inv_w2 = inv_w2;

end