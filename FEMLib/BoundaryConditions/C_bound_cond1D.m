function [A,b,u_g] = C_bound_cond1D(A,b,femregion,Dati)
%% [A,b,u_g] = C_bound_cond1D(A,b,femregion,Dati)
%==========================================================================
% Assign Dirchlet boundary conditions
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          femregion   : (struct)  see C_create_femregion.m
%          Dati        : (struct)  see C_dati.m

%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          u_g         : (sparse(ndof,1) real) evaluation of Dirichlet conditions 
%

%fprintf('============================================================\n')
%fprintf('Assign Dirichlet boundary conditions ... \n');
%fprintf('============================================================\n')

ndof = length(b);
u_g = sparse(ndof,1);

if strcmp(Dati.bc1, "D") && strcmp(Dati.bc2, "D")
    boundary_points = femregion.boundary_points;
elseif strcmp(Dati.bc1, "D")
    boundary_points = femregion.boundary_points(1);
else
    boundary_points = femregion.boundary_points(2);
end

% x = femregion.dof(boundary_points,1);
% t = Dati.t;
u_g(boundary_points) = 0; % Compute the lifting operator ug

x_g = sparse(ndof,1);
A_0 = A;

b_0 = b-A*u_g; % modify the load vector --> F(v) = F(v) - a(ug,v)

% Reduce the system A in order to solve the pb with 
% homogeneous Dirichlet conditions 
for k = 1:length(boundary_points)
    A_0(boundary_points(k),:) = 0;
    A_0(:,boundary_points(k)) = 0;
    A_0(boundary_points(k),boundary_points(k)) = 1;
    b_0(boundary_points(k)) = 0;
end

b = b_0;
A = A_0;
