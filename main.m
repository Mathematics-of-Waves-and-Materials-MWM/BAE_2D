%%%%% load the geometry of the obstacle
input_geometry

 %body_boundary = square; %%%% a square obstacle
 % N_angles = 8; %%% number of independent solutions for embedding

% body_boundary = strip; %%%% a square obstacle
% N_angles = 4; %%% number of independent solutions for embedding

% body_boundary = P_figure; %%%% a square obstacle
% N_angles = 10; %%% number of independent solutions for embedding

  body_boundary = g_strip; %%%% a square obstacle
 N_angles = 6; %%% number of independent solutions for embedding

 %%%%%%% build a mesh; the obstacle should belong to domain [-50,50]x[-50,50]

 mesher


 %%%% Assemble FEM matrices K and M 

 Assemble_K_M

 %%% solve boundary algebraic equations

 solver

 %%% recover the field everywhere

 field_recovery_for_plotting

 %%% calculate direcitivities, check embedding formula

 embedding



