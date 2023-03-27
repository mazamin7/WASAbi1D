function C = get_residue_matrix(N_x)
    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    % delta = -49/18; % unused

    C = sparse(N_x,N_x);

    C(N_x/2-2,N_x/2:N_x/2+1) = [-alpha, alpha];
    C(N_x/2-1,N_x/2-1:N_x/2+2) = [-alpha, -beta, beta, alpha];
    C(N_x/2,N_x/2-2:N_x/2+3) = [-alpha, -beta, -gamma, gamma, beta, alpha];

    C(N_x/2+1,N_x/2-2:N_x/2+3) = -[-alpha, -beta, -gamma, gamma, beta, alpha];
    C(N_x/2+2,N_x/2-1:N_x/2+2) = -[-alpha, -beta, beta, alpha];
    C(N_x/2+3,N_x/2:N_x/2+1) = -[-alpha, alpha];
end
