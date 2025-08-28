if sum(mask_squares)<0.5

%   connectivity_adjacent =  connectivity(mask_adjacent_squares,:);


%%%%%%% group the adjacent connectivities into two sets laying on different
%%%%%%% sides of the boundary
   prev_index = [];
   cur_indices = [];
   mask_con_up = false(num_squares,1);
   cur_bound_connectivity = zeros(1,2);
   for n_cur = 1: size(body_boundary,1)-1
       cur_mask_squares = true(num_squares,1);
       for m_cur = 0:1
        idx = boundary_idx_ar(n_cur+m_cur);  
        cur_bound_connectivity(m_cur+1) = idx;
        connectivity_dist = abs(connectivity - idx);
        connectivity_dist = min(connectivity_dist,[],2);
        cur_mask_squares = cur_mask_squares&connectivity_dist<0.5;
        cur_indices = find(cur_mask_squares);
       end
       if isempty(prev_index)
          mask_con_up(cur_indices(2)) = true;
          prev_index = cur_indices(2);
       else
           prev_con = connectivity(prev_index,:);
           cur_con  = connectivity(cur_indices(1),:);
           prev_con = sort(prev_con);
           cur_con = sort(cur_con);
           delta = intersect(prev_con,cur_con);
           flag = length(delta);
           if flag>1.5
              mask_con_up(cur_indices(1)) = true;
              prev_index = cur_indices(1);
           else
              mask_con_up(cur_indices(2)) = true;
              prev_index = cur_indices(2);
           end
       end
   
   end

%%%%% duplicate nodes

num_add_nodes = length(boundary_idx_ar(2:end-1));

node_coords = [node_coords;node_coords(boundary_idx_ar(2:end-1),:)];

mask_boundary = [mask_boundary, true(1,num_add_nodes)];
mask_outer_nodes = [mask_outer_nodes; true(num_add_nodes,1)];
mask_adjacent_boundary = [mask_adjacent_boundary;num_nodes + 1:num_add_nodes];


con_indices = find(mask_con_up);
%%%%% adjust one set of adjacent connectivities
   for n_cur  = 1 : sum(mask_con_up)
       cur_connectivity = connectivity(con_indices(n_cur),:);
       [nodes_to_change,ia,ib] = intersect(cur_connectivity,boundary_idx_ar(2:end-1));
       cur_connectivity(ia) = ib + num_nodes;
       connectivity(con_indices(n_cur),:) = cur_connectivity;
   end

   num_nodes = size(node_coords,1);

   % figure;
   %  patch('Faces', connectivity(mask_con_up,:), 'Vertices', node_coords, ...
   %     'FaceColor', 'none', ...        % No face fill
   %     'EdgeColor', 'k');              % Black edges
   % axis equal
   %    hold all
   % plot(node_coords(mask_boundary,1),node_coords(mask_boundary,2),'r*')
end



