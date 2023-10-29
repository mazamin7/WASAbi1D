function K = get_stiffness_matrix(N_x, space_order)

    assert(N_x == 2*floor(N_x/2), 'N_x is odd');

    assert(space_order == 2 || space_order == 4 || space_order == 6 || space_order == 8);
    n = space_order/2;

    if n == 1
        coefs = [1, -2];
    elseif n == 2
        coefs = [-1/12, 4/3, -5/2];
    elseif n == 3
        coefs = [1/90, -3/20, 3/2, -49/18];
    else
        coefs = [-1/560, 8/315, -1/5, 8/5, -205/72];
    end

    K = sparse(N_x,N_x);

    for i = 1:n
        row = i;

        col1 = 1;
        col2 = row+n;

        vals = [coefs(n+2-i:n+1) coefs(n:-1:1)];
        K(row,col1:col2) = vals;

        row = N_x-i+1;

        col1 = row-n;
        col2 = N_x;

        vals = [coefs(1:n) coefs(n+1:-1:n+2-i)];
        K(row,col1:col2) = vals;
    end

    for i = n+1:N_x-n
        row = i;
        
        col1 = i-n;
        col2 = i+n;

        vals = [coefs(1:n+1) coefs(n:-1:1)];
        K(row,col1:col2) = vals;
    end

    C = get_residue_matrix(N_x, space_order);

    K(1:n,1:n) = K(1:n,1:n) - C(N_x/2+1:N_x/2+n,N_x/2+1:N_x/2+n);
    K(N_x-n+1:N_x,N_x-n+1:N_x) = K(N_x-n+1:N_x,N_x-n+1:N_x) - C(N_x/2+1-n:N_x/2,N_x/2+1-n:N_x/2);

end
