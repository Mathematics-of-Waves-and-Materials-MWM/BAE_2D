% clear all
% load('mesh.mat')

%%%% element matrices for finite difference method

K_el = -1/2*[2,-1,0,-1;-1,2,-1,0;0,-1,2,-1;-1,0,-1,2];
M_el = eye(4)/4;
Embedding_el = [0,0.5,0,0;0.5,0,0,0;0,0,0,0.5;0,0,0.5,0];


Kmat = sparse(num_nodes,num_nodes);
Mmat = sparse(num_nodes,num_nodes);
Emat = sparse(num_nodes,num_nodes);

for n_el = 1:size(connectivity,1)
    n_el
    
    nodes = connectivity(n_el , :) ;  % 1 x 4 


    
    %%%%%%%%%%%%%%%%%%%%%%%%%

    
         Kmat(nodes,nodes) = Kmat(nodes,nodes) + K_el ; 
         Mmat(nodes,nodes) = Mmat(nodes,nodes) + M_el ;
         Emat(nodes,nodes) = Emat(nodes,nodes)+Embedding_el;
         
end 

Emat(Emat>0) = 1;
Kbound = Kmat(mask_boundary,mask_adjacent_boundary);
Mbound = Mmat(mask_boundary,mask_adjacent_boundary);

Kfull = Kmat(mask_boundary,mask_outer_nodes);
Mfull = Kmat(mask_boundary,mask_outer_nodes);


save('FEM_matrices','Kbound','Mbound','Kfull','Mfull','Kmat','Mmat','Emat')