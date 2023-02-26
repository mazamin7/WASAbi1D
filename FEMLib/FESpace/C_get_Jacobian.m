function [BJ, pphys_1D] = C_get_Jacobian(loc_coord, nodes_1D)
%% [BJ, pphys_1D] = C_get_Jacobian(loc_coord, nodes_2D)
%==========================================================================
% Compute Jacobian of the map transformation and physical nodes
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          loc_coord   : (array real) coordinates of grid nodes
%          nodes_1D    : (struct) see C_quadrature.m
%
%    OUTPUT:
%          BJ          :  (array real) [2 x 2 x nNode] Jacobian
%          pphys_1D    :  (array real) map of the reference element nodes

for k = 1:length(nodes_1D)
    
    x0 = loc_coord(1,1);   % x-coordinates of vertices
    x1 = loc_coord(2,1);
 
    BJ = x1-x0;   % Jacobian of elemental map
       
    pphys_1D(k) = x0 + BJ * nodes_1D(k);
        
end


