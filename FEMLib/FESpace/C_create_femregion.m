function [femregion] = C_create_femregion(Dati,Region) 
%% [femregion] = C_create_femregion(Dati,Region)
%==========================================================================
% Creates conforming finite element space
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          Region      : (struct)  see C_create_mesh.m
%
%    OUTPUT:
%          femregion    : (struct) 

degree = 1;
nln = degree + 1;  
        
[bound_pts] = C_create_bound_pts(Dati.domain,Region.coord);

%==========================================================================
% COSTRUZIONE STRUTTURA FEMREGION
%==========================================================================
femregion=struct('fem',Dati.fem,...
                'domain',Region.domain,...
                'h',Region.h,...
                'nln',nln,...
                'ndof',length(Region.coord),...
                'ne',Region.ne,...
                'dof',Region.coord,...
                'nqn_1D',Dati.nqn_1D,...
                'degree',degree,...
                'coord',Region.coord,...
                'connectivity',Region.connectivity,...
                'boundary_points',bound_pts);         
            
function [bound_pts] = C_create_bound_pts(domain, coord)
%% [bounds_pts] = C_create_bound_pts(domain, coord)
%==========================================================================
% Creates boundary point list
%==========================================================================
%    called in C_create_femregion.m
%
%    INPUT:
%          domain     : [2 x 1] (array real)  extrema of the rectangle
%          coord      : [Np x 1] (array real) coordinates of the grid
%                                             points
%
%    OUTPUT:
%          bound_pts  : [Nbp x 1] (array int) indices of the boundary
%                                             points

x0 = domain(1,1);
x1 = domain(1,2);
bound_pts = ones(length(coord),1);

bound_pts = find(coord(:,1)==x0 | coord(:,1)==x1 );

        
            
            
