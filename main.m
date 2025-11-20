clear all
%%%%%%%%% Launches all the scripts alltogether
%%%%% Note that some of the scripts could be launched only once

%%%%% load the geometry of the obstacle
input_geometry

%%%%% sample geometries to choose

 geometry = geometry_square; %%%% a square obstacle
 N_angles = 8; %%% number of independent solutions for embedding

% geometry = geometry_strip; %%%% a strip
% N_angles = 4; %%% number of independent solutions for embedding

% geometry = geometry_P_figure; %%%% a P shaped figure
% N_angles = 10; %%% number of independent solutions for embedding

% geometry = geometry_g_strip; %%%% g-shaped thin line
% N_angles = 6; %%% number of independent solutions for embedding

% geometry = geometry_square_grid; %%%% a grid of 9 squares
% N_angles = 9*8; %%% number of independent solutions for embedding

% geometry = geometry_open_square; %%%% a grid of 9 squares
% N_angles = 10; %%% number of independent solutions for embedding



 %%%%%%% build a mesh; the obstacle should belong to domain [-50,50]x[-50,50]

 mesher


 %%%% Assemble FEM matrices K and M 

 Assemble_K_M

 %%%%%%%%% set wavenumber
 K = 0.3 + 0.01i;
 %%%%%%% set incident parameter (cotanget of the angle of incidence)


 %%%%%% Compute free field Greens function
 Compute_Greens_function

 %%% solve boundary algebraic equations

 beta_star = 0.9;

 solver

 %%% recover the field everywhere

 field_recovery_for_plotting

 %%% calculate direcitivities, check embedding formula

 embedding



