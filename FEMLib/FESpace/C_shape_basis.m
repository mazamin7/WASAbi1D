function [basis] = C_shape_basis(Dati)
%% function [basis] = C_shape_basis(Dati)
%==========================================================================
% Compute the shape functions for triangular and quadrilateral elements
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%
%    OUTPUT:
%          basis       : (struct) .num  (int) number of basis functions
%                                 .n_edge (int) number of end-points
%                                 .fbases (string) basis functions
%                                 .Gbases (string) df/dx

% fprintf('============================================================\n')
% fprintf('Compute %s finite elements...\n',Dati.fem);
% fprintf('============================================================\n')

nln = 2;
basis = struct('num',nln,...
               'n_edge',2,...
               'fbases',{'1-csi',...
                         'csi',...
                        },...
               'Gbases',{'-1+0.*csi',...
                         ' 1+0.*csi',...
                        });
        
        
