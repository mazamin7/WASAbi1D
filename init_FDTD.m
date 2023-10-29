function data = init_FDTD(len_x, c, dt, dh, space_order, temp_order, alpha_abs)

    N = floor(len_x/dh);

    A = get_stiffness_matrix(N, space_order);

    data.N = N;
    data.A = A;
    data.c = c;
    data.dt = dt;
    data.dh = dh;
    data.alpha_abs = alpha_abs;
    data.temp_order = temp_order;

end