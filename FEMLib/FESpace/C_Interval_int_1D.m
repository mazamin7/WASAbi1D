function [node_1D,w_1D]=C_Interval_int_1D(nqn_1D)
%% [node_1D,w_1D]=C_Interfval_int_1D(nqn_1D)
%==========================================================================
% Compute Gauss-Legendre nodes and weights on the reference triangle
%==========================================================================
%    called in C_quadrature.m
%
%    INPUT:
%          nqn_1D      : (integer) number of quadrature nodes
%
%    OUTPUT:
%          node_1D     : (array) [nqn_2D x 1] list of quadrature nodes
%          w_1D        : (array) [nqn_2D x 1] list of quadrature weights        
       
switch nqn_1D 
    
case{1}
    %===================================================
    % mid point rule
    %===================================================
    
    xnod = 0.5;
    w = 1;
    
case{2}
    
    %===================================================
    % Trapezi rule
    %===================================================
    
    w(1) = 0.5;
    w(2) = 0.5;
    
    xnod(1) = 0;
    xnod(2) = 1;

    
case{3}
    %===================================================
    % Cavalieri-Simpson rule
    %===================================================
        
    w(1) = 1/3;
    w(2) = 4/3;
    w(3) = 1/3;

    xnod(1) = 0;
    xnod(2) = 0.5;
    xnod(3) = 1;
    
end



node_1D=[xnod'];
w_1D=w;    
