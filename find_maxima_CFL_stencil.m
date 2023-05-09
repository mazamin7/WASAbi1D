function max_amplitude = find_maxima_CFL_stencil()
    alpha = 1/90;
    beta = -3/20;
    gamma = 3/2;

    % Define the function
    f = @(x) alpha*sin(3*x).^2 + beta*sin(2*x).^2 + gamma*sin(x).^2;
    
    % Find the maximum amplitude using fminbnd
    [~, max_amplitude] = fminbnd(@(x) -f(x), 0, pi);
    
    % Return the absolute value of the maximum amplitude
    max_amplitude = abs(max_amplitude);
end
