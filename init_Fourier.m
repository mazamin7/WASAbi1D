function Fourier_data = init_Fourier(len_x, c, dt, dh, order, alpha)

    % modal analysis
    inv_lx = 1 / len_x;

    N = floor(len_x/dh);
    
    cwt = zeros(N, 1);
    w2 = zeros(N, 1);
    swt = zeros(N, 1);
    w = zeros(N, 1);
    inv_w = zeros(N, 1);
    inv_w2 = zeros(N, 1);

    alpha_abs = alpha;
    alpha2 = alpha_abs * alpha_abs;
    eatm = exp(- alpha_abs * dt);
    
	% first element is empty but no one cares
	n = 2:N;
	
    w_over_pi = (n-1) * c * inv_lx;
    w(n) = pi * w_over_pi;
    inv_w(n) = 1./w(n);
    inv_w2(n) = inv_w(n) .* inv_w(n);
    w2(n) = w(n) .* w(n);
    cwt(n) = cospi(w_over_pi * dt);
    swt(n) = sinpi(w_over_pi * dt);

    Fourier_data.N = N;
    Fourier_data.w2 = w2;
    Fourier_data.cwt = cwt;
    Fourier_data.dt = dt;
    Fourier_data.order = order;
    Fourier_data.swt = swt;
    Fourier_data.alpha_abs = alpha_abs;
    Fourier_data.alpha2 = alpha2;
    Fourier_data.eatm = eatm;
    Fourier_data.w = w;
    Fourier_data.inv_w = inv_w;
    Fourier_data.inv_w2 = inv_w2;

end