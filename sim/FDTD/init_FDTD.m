function data = init_FDTD(len_x, c, dt, dh, space_order, temp_order, alpha_abs)

    laplacian = get_laplacian_kernel(space_order);

    data.laplacian = laplacian;
    data.c = c;
    data.dt = dt;
    data.dh = dh;
    data.alpha_abs = alpha_abs;
    data.temp_order = temp_order;

end