function FEM_data = init_FEM(N,c,dt,dh,isDamped,alpha_abs,boundCond1,boundCond2)
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
    addpath FEMLib/Assembly
    addpath FEMLib/BoundaryConditions
    addpath FEMLib/MeshGeneration
    addpath FEMLib/FESpace
    
    %==========================================================================
    % LOAD DATA FOR TEST CASE
    %==========================================================================
    
    Dati = C_dati();

    Dati.N = N;

    Dati.domain(2) = (N - 1) * dh;

    c2 = c^2;
    Dati.c2 = c2;
    Dati.dt = dt;

    Dati.bc1 = boundCond1;
    Dati.bc2 = boundCond2;
    
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

    if(strcmp(Dati.bc1,'R'))
        A_nbc(1,1)     = A_nbc(1,1) + Dati.c2;
    end

    if(strcmp(Dati.bc2,'R'))
        A_nbc(end,end) = A_nbc(end,end) + Dati.c2;
    end
    
    FEM_data.M_nbc = M_nbc;
    FEM_data.A_nbc = A_nbc;
    FEM_data.alpha_abs = alpha_abs;
    FEM_data.isDamped = isDamped;
    FEM_data.Dati = Dati;
    FEM_data.femregion = femregion;
    
end


