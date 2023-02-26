%=======================================================================================================
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================

function [Dati]=C_dati()

Dati = struct( 'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'N', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing          
               'fem',              'P1',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS' ...        
               ... % uniform regular mesh
               );



