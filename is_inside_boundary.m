function state= is_inside_boundary(node,boundary)

% determine whether the current node is strictly inside the boundary

   state_right = check_to_the_right(node,boundary);
   state_left =  check_to_the_right(node([2,1]),boundary(:,[2,1]));

   state = state_right&state_left;


end

function   state = check_to_the_right(node,boundary) 
       max_node_m = max(boundary(:,1));
    
       mask_boundary = boundary == node;
       mask_boundary = mask_boundary(:,1)&mask_boundary(:,2);
    
       if(sum(mask_boundary)>0.5)
           state = false(1,1);
       else 
           half_line_m = node(1):max_node_m;
           N_nodes = length(half_line_m);
           half_line = [half_line_m; repmat(node(2),1,N_nodes)].';
           intersection_index = 0;
           for n_cur = 1:N_nodes
               cur_node = half_line(n_cur,:);
               mask_boundary = boundary == cur_node;
               mask_boundary = mask_boundary(:,1)&mask_boundary(:,2);
               if (sum(mask_boundary)>0.5)
                  intersection_index = intersection_index + 1;
               end
           end
    
           if(mod(intersection_index,2)>0.5)
           state = true(1,1);
           else 
            state = false(1,1);   
           end
       end
end