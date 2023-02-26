% quadrature nodes and whieghs for faces integrals
function  [node_1D,w_1D] = C_quadrature(Dati)
%% function [node_1D,w_1D] = C_quadrature(Dati)
%==========================================================================
% Compute quadrature nodes and weights triangular and quadrilateral elements
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%
%    OUTPUT:
%          node2D      : (struct) .num  (int) number of basis functions
%                                 .nedge (int) number of edges
%                                 .fbases (string) basis functions
%                                 .Gbases  (string) df/dx
%

[node_1D,w_1D] = C_Interval_int_1D(Dati.nqn_1D);
        
