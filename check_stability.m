function [stable] = check_stability(c, dt, dh, alpha_abs, order, xi, nu)

    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;
    delta = -49/18;


    % check stability
    if order == 1
        % expected max eig of K
        max_eig_K_abs = abs(alpha)+abs(beta)+abs(gamma)+abs(delta)+abs(gamma)+abs(beta)+abs(alpha);
        
        
        % expected max eig of Koverline
        max_eig_Koverline1 = -alpha_abs + 1i*sqrt(-alpha_abs^2 + c^2/dh^2*max_eig_K_abs);
        max_eig_Koverline2 = -2*alpha_abs;
        lim = 1/2*sqrt(c^2/dh^2*max_eig_K_abs);
        
        if alpha_abs < lim
            max_eig_Koverline = max_eig_Koverline1;
        else
            max_eig_Koverline = max_eig_Koverline2;
        end
        
        
        % expected max eig of B
        max_eig_B1 = -dt * nu * max_eig_Koverline + sqrt((dt * nu * max_eig_Koverline)^2 + xi);
        max_eig_B2 = -dt * nu * max_eig_Koverline - sqrt((dt * nu * max_eig_Koverline)^2 + xi);
        max_eig_B3 = xi;
        
        [max_eig_B, ~] = complex_max(max_eig_B1, max_eig_B2, max_eig_B3);
    
        stable = abs(max_eig_B) < 1; % check that largest abs is less than 1
    else
        var_1511 = 1/2 * (abs(alpha) + abs(beta) + abs(gamma) - delta/2);
        CFL = 1 / sqrt(var_1511);
        stable = dt < dh / c * CFL; % CFL condition
    end

end