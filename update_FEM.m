function p_next = update_FEM(p_curr,p_prev,c,dt,dh,force)
%==========================================================================
% Solution of the Wave Equation with linear finite elements 
% coupled with the leap-frog scheme 
% 
%  u_tt - c^2 u_xx = f  in (a,b)x(0,T)
%  u_t(x,0) = v_0       in (a,b) 
%  u(x,0)   = u_0       in (a,b)
%  + boundary conditions:
%  Dirichlet:  u(s,t) = g  s=a,b
%  Neumann  : c^2du/dn(s,t) = g s=a,b
%  Periodic : u(a,t) = a(b,t)
%  Absorbing: du/dt(s,t) + cdu/dn(s,t) = 0  s=a,b 
%==========================================================================
%
%    INPUT:
%          TestName    : (string)  see C_dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Dati        : (struct)  see C_dati.m
%
% Usage:
%    [femregion,Dati] = C_main1D('Test1')

    addpath FEMLib
%     addpath Assembly
%     addpath BoundaryConditions
%     addpath MeshGeneration
%     addpath FESpace
    
    %==========================================================================
    % LOAD DATA FOR TEST CASE
    %==========================================================================
    
    Dati = C_dati();

    N = length(p_curr);
    Dati.N = N;

    Dati.domain(2) = N * dh;

    c2 = c^2;
    Dati.c2 = c2;
    Dati.dt = dt;
    Dati.force = force;
    
    %==========================================================================
    % MESH GENERATION
    %==========================================================================
    
    [Region] = C_create_mesh(Dati);
    
    %==========================================================================
    % FINITE ELEMENT REGION
    %==========================================================================
    
    [femregion] = C_create_femregion(Dati,Region);
    
    %==========================================================================
    % BUILD FINITE ELEMENT MATRICES
    %==========================================================================
    
    [M_nbc,A_nbc] = C_matrix1D(Dati,femregion);

    if(strcmp(Dati.bc,'R'))
        A_nbc(1,1)     = A_nbc(1,1) + Dati.c2;
        A_nbc(end,end) = A_nbc(end,end) + Dati.c2;
    end
    
    u1 = p_curr;
    u0 = p_prev;

    %==========================================================================
    % BUILD FINITE ELEMENTS RHS a time t
    %==========================================================================
    [b_nbc] = C_rhs1D(Dati,femregion);
    
    if(strcmp(Dati.bc,'N'))
        b_nbc(1)   = b_nbc(1);%   + eval(Dati.neumann1);
        b_nbc(end) = b_nbc(end);% + eval(Dati.neumann2);
    elseif(strcmp(Dati.bc,'R'))
        b_nbc(1)   = b_nbc(1);%   - Dati.c2*eval(Dati.neumann1);
        b_nbc(end) = b_nbc(end);% + Dati.c2*eval(Dati.neumann2);
    end
    
    % Repeat steps 1) to 5) for the general time step
    b_nbc = Dati.dt^2 * (b_nbc - A_nbc*u1) + 2 * M_nbc * u1 - M_nbc * u0;
    
    if(strcmp(Dati.bc,'D'))
        [M,b,u_g] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
        u2 =  M\b;
        u2 = u2 + u_g;
    elseif(strcmp(Dati.bc,'N') || strcmp(Dati.bc,'R'))
        u2 =  M_nbc\b_nbc;
    elseif(strcmp(Dati.bc,'P'))
        b_nbc(1,:) = b_nbc(1,:) + b_nbc(end,:);
        b_nbc(end) = 0;
        u2 =  M_nbc\b_nbc;
    elseif(strcmp(Dati.bc,'A'))
        b_nbc(1)   = b_nbc(1)   - sqrt(Dati.c2)*(u1(1)-u0(1))*Dati.dt;
        b_nbc(end) = b_nbc(end) - sqrt(Dati.c2)*(u1(end)-u0(end))*Dati.dt;
        u2 =  M_nbc\b_nbc;
    end
    
    p_next = u2;
    
end


