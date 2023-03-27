function FEM_data = init_FEM(len_x,c,dt,dh,isDamped,alpha_abs,boundCond1,boundCond2)

    % Loading FEMLib dependencies
    addpath FEMLib
    addpath FEMLib/Assembly
    addpath FEMLib/BoundaryConditions
    addpath FEMLib/MeshGeneration
    addpath FEMLib/FESpace
    
    %==========================================================================
    % LOAD DATA
    %==========================================================================
    
    Dati = C_dati();

    N = floor(len_x/dh);
    N = 2 * floor(N/2);
    Dati.N = N;

    Dati.domain(2) = len_x;

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

    % Storing data in data structure
    FEM_data.M_nbc = M_nbc;
    FEM_data.A_nbc = A_nbc;
    FEM_data.alpha_abs = alpha_abs;
    FEM_data.isDamped = isDamped;
    FEM_data.Dati = Dati;
    FEM_data.femregion = femregion;
    
end


