% clear all
% 
% load('mesh.mat')
% load('FEM_matrices.mat')
%load('Green_function.mat')




%%% Assemble propagator matrix

Amat = zeros(num_nodes_boundary);
Amat_adjacent = zeros(num_nodes_adjacent,num_nodes_boundary);
Afull = zeros(num_nodes_outer,num_nodes_boundary);

for n_cur = 1:num_nodes_boundary

    cur_nodes = abs(node_coords_boundary(:,:) - node_coords_boundary(n_cur,:)) + ones(num_nodes_boundary,2);
    linear_idx = sub2ind(size(Greens_func), cur_nodes(:, 1), cur_nodes(:, 2));
    Amat(n_cur,:) = Greens_func(linear_idx);
   
    cur_nodes = abs(node_coords_adjacent(:,:) - node_coords_boundary(n_cur,:)) + ones(num_nodes_adjacent,2);
    linear_idx = sub2ind(size(Greens_func), cur_nodes(:, 1), cur_nodes(:, 2));
    Amat_adjacent(:,n_cur) = Greens_func(linear_idx);

    cur_nodes = abs(node_coords_outer(:,:) - node_coords_boundary(n_cur,:)) + ones(num_nodes_outer,2);
    linear_idx = sub2ind(size(Greens_func), cur_nodes(:, 1), cur_nodes(:, 2));
    Afull(:,n_cur) = Greens_func(linear_idx);

    
end

 %%% recover the field everywhere 

Afull_adj = zeros(num_nodes_adjacent,num_nodes_outer);

 for n_cur = 1: num_nodes_adjacent
    cur_nodes = abs(node_coords_outer(:,:) - node_coords_adjacent(n_cur,:)) + ones(num_nodes_outer,2);
    linear_idx = sub2ind(size(Greens_func), cur_nodes(:, 1), cur_nodes(:, 2));
    Afull_adj(n_cur,:) = Greens_func(linear_idx);
 end



%%% Assemble rhs matrices


cur_nodes = abs(node_coords_boundary(:,:)) + ones(num_nodes_boundary,2);
linear_idx = sub2ind(size(Amat), cur_nodes(:, 1), cur_nodes(:, 1));
G_rhs = Greens_func(linear_idx);

%% incidence parameters beta s_in q_in
%N_angles = 8;
dtheta_ar = pi/2/N_angles; 
theta_ar_in = 0.1:dtheta_ar:pi/2;

beta_ar = tan(theta_ar_in);

sol_ar = zeros(num_nodes_boundary,N_angles);
s_ar = zeros(1,N_angles);
q_ar = zeros(1,N_angles);
for n_cur = 1:N_angles

[s1_in,s2_in,s3_in,s4_in] = saddle_points_torus(K,beta_ar(n_cur));

s_ar_cur = [s1_in,s2_in,s3_in,s4_in];

s_in = s_ar_cur(abs(s_ar_cur)>1&abs(s_ar_cur)<1.2);

q_in = -(K^2 - 4 + s_in + 1./s_in)/2 + sqrt((K^2-4+s_in + 1./s_in).^2 - 4)/2;

if abs(q_in)<1
  q_in = 1/q_in;
end

%(s_in - 1./s_in)/(q_in-1./q_in);

s_ar(n_cur) = s_in;
q_ar(n_cur) = q_in;

uin = s_in.^(node_coords_boundary(:,1)).*q_in.^(node_coords_boundary(:,2));

F_ar = - Amat_adjacent.'*(Kbound + K^2*Mbound).'*uin;

sol_ar(:,n_cur) = Amat\F_ar;

end

%%%%% one more computation for target value
%%%%% incidence parameter to calculate
 theta_star = atan(beta_star);
 [s_star,q_star] = find_plane_wave_parameters(K,beta_star);
uin = s_star.^(node_coords_boundary(:,1)).*q_star.^(node_coords_boundary(:,2));

F_ar = - Amat_adjacent.'*(Kbound + K^2*Mbound).'*uin;

sol_star = Amat\F_ar;


 save('result','sol_ar','K','beta_ar','Afull_adj','Afull','q_ar','s_ar',...
     'theta_ar_in','N_angles','s_star','q_star','beta_star','theta_star','sol_star')











