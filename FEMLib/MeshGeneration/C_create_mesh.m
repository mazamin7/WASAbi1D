function [Region] = C_create_mesh(Dati)
%% [Region] = C_create_mesh(Dati)
%==========================================================================
% Creates triangular or quadrilateral mesh
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%
%    OUTPUT:
%          Region      : (struct) having fields: dimension
%                                                domain 
%                                                mesh size
%                                                number of vertices
%                                                number of elements
%                                                coordinates
%                                                boundary points
%                                                connectivity

%

x0 = Dati.domain(1);
xL = Dati.domain(2);

%================================================
% GEOMETRICAL INFO
 nEl = Dati.N - 1; 
 nVert = nEl + 1;
 p = linspace(x0,xL,nVert);
 t = [[1:nVert-1]' [2:nVert]']';
 MeshSize = (xL-x0)./nEl;
%================================================

% struttura dati della mesh
Region = struct('dim',1,...
               'domain',Dati.domain,...
               'h',MeshSize,...
               'nvert',nVert,...
               'ne',nEl,...
               'coord',p',...
               'boundary_points',[x0,xL],...
               'connectivity',t);
           
           
