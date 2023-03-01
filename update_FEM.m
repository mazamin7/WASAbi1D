function p_next = update_FEM(FEM_data,p_curr,p_prev,force)
% Computes p_next given p_curr, p_prev, force and FEM_data
%
% Inputs:
%   - FEM_data: data structure containing data needed for the
%   computation
%   - p_curr: the current pressure values (a column vector)
%   - p_prev: the previous pressure values (a column vector)
%   - force: the applied force (a column vector)
%
% Output:
%   - p_next: the next pressure values (a column vector)

    addpath FEMLib
    addpath FEMLib/Assembly
    addpath FEMLib/BoundaryConditions
    addpath FEMLib/MeshGeneration
    addpath FEMLib/FESpace
    
    %==========================================================================
    % LOAD DATA
    %==========================================================================
    
    M_nbc = FEM_data.M_nbc;
    A_nbc = FEM_data.A_nbc;
    alpha_abs = FEM_data.alpha_abs;
    isDamped = FEM_data.isDamped;
    Dati = FEM_data.Dati;
    femregion = FEM_data.femregion;

    Dati.force = force;
    
    u1 = p_curr;
    u0 = p_prev;

    %==========================================================================
    % BUILD FINITE ELEMENTS RHS at current time instant
    %==========================================================================
    [b_nbc] = C_rhs1D(Dati,femregion);
    
    if(strcmp(Dati.bc1,'N'))
        b_nbc(1)   = b_nbc(1);%   + eval(Dati.neumann1);
    elseif(strcmp(Dati.bc2,'R'))
        b_nbc(1)   = b_nbc(1);%   - Dati.c2*eval(Dati.neumann1);
    end

    if(strcmp(Dati.bc2,'N'))
        b_nbc(end) = b_nbc(end);% + eval(Dati.neumann2);
    elseif(strcmp(Dati.bc2,'R'))
        b_nbc(end) = b_nbc(end);% + Dati.c2*eval(Dati.neumann2);
    end
    
    % Repeat steps 1) to 5) for the general time step
    if isDamped == false
        b_nbc = Dati.dt^2 * (b_nbc - A_nbc*u1) + 2 * M_nbc * u1 - M_nbc * u0;
    else
        b_nbc = Dati.dt^2 * (b_nbc - A_nbc*u1) ...
            - Dati.dt * alpha_abs * M_nbc * (u1 - u0) ...
            + 2 * M_nbc * u1 - M_nbc * u0;
    end
    
    if(strcmp(Dati.bc1,'D') || strcmp(Dati.bc2,'D'))
        [M,b,u_g] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
        u2 =  M\b;
        u2 = u2 + u_g;
    else
        if strcmp(Dati.bc1,'A')
            b_nbc(1)   = b_nbc(1)   - sqrt(Dati.c2)*(u1(1)-u0(1))*Dati.dt;
        end

        if strcmp(Dati.bc2,'A')
            b_nbc(end) = b_nbc(end) - sqrt(Dati.c2)*(u1(end)-u0(end))*Dati.dt;
        end

        u2 =  M_nbc\b_nbc;
    end
    
    p_next = u2;
    
end


