% clear all
% 
% load('mesh.mat')
% load("result.mat")
% load('FEM_matrices.mat')

n_cur = 2;
sol = sol_ar(:,n_cur);
s_in = s_ar(n_cur);
q_in = q_ar(n_cur);
uin = s_in.^(node_coords_boundary(:,1)).*q_in.^(node_coords_boundary(:,2));

%% reassemble field for illustrations

scat_field = Afull*sol + Afull_adj.'*(Kbound + K^2*Mbound).'*uin;


pressure = zeros(num_nodes,1);

pressure(mask_outer_nodes) = scat_field;
pressure(mask_boundary) = -uin;

uin_everywhere = s_in.^(node_coords(:,1)).*q_in.^(node_coords(:,2));

pressure = pressure + uin_everywhere;

fig = figure;
patch('Faces', connectivity, 'Vertices', node_coords, 'FaceVertexCData', real(pressure), ...
      'FaceColor', 'interp');
axis equal
   hold all
   plot(node_coords(boundary_idx_ar,1),node_coords(boundary_idx_ar,2),'k-','LineWidth',2)
 shading interp
 title('${\rm Re}[u]$',FontSize=18,Interpreter='latex')
 xlabel('$m$',FontSize=16,Interpreter='latex')
 ylabel('$n$',FontSize=16,Interpreter='latex')
 colorbar

 % exportgraphics(fig, 'field.pdf', ...
 %    'ContentType', 'vector')

Kmat(mask_boundary,:)=0;
Mmat(mask_boundary,:)=0;
Kmat(mask_boundary,mask_boundary)=eye(num_nodes_boundary);

delta = (Kmat + K^2*Mmat)*pressure;

%%%% delta = (Kmat + K^2*Mmat)*uin_everywhere;


 % figure;
 % patch('Faces', connectivity, 'Vertices', node_coords, 'FaceVertexCData', abs(delta), ...
 %       'FaceColor', 'interp');
 % axis equal
 %    hold all
 %    plot(node_coords(mask_boundary,1),node_coords(mask_boundary,2),'r*')

% figure; plot(sort(abs(delta(mask_outer_nodes))))

%

%% embedding games
% 
 Hp = Emat*pressure - (s_in + 1/s_in)*pressure;
% 
 delta = (Kmat + K^2*Mmat)*Hp;


% figure;
% patch('Faces', connectivity, 'Vertices', node_coords, 'FaceVertexCData', abs(delta), ...
%       'FaceColor', 'interp');
% axis equal
%    hold all
%    plot(node_coords(mask_boundary,1),node_coords(mask_boundary,2),'r*')