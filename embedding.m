% clear all
% 
% load('mesh.mat')
% load("result.mat")
% load('FEM_matrices.mat')



%%%% This script calculates directivity using embedding formula and
%%%% compares the result with direct computation
%% Compute directivity

theta_ar_obs = 0.015:0.017:pi/2-0.015;
N_angles_obs = length(theta_ar_obs);
beta_ar_obs = tan(theta_ar_obs);

s_ar_obs = zeros(1,N_angles_obs);
q_ar_obs = zeros(1,N_angles_obs);

for n_cur = 1:N_angles_obs 

[s1_in,s2_in,s3_in,s4_in] = saddle_points_torus(K,beta_ar_obs(n_cur));

s_ar_cur = [s1_in,s2_in,s3_in,s4_in];

s_in = s_ar_cur(abs(s_ar_cur)>1&abs(s_ar_cur)<1.2);

if(length(s_in)>1)
    [~,ind] = max(abs(imag(s_in)));
    s_in = s_in(ind);
end



q_in = -(K^2 - 4 + s_in + 1./s_in)/2 + sqrt((K^2-4+s_in + 1./s_in).^2 - 4)/2;

if abs(q_in)<1
  q_in = 1/q_in;
end

s_ar_obs(n_cur) = s_in;
q_ar_obs(n_cur) = q_in;

end


 

   phase = s_ar_obs.^(node_coords_boundary(:,1)).* q_ar_obs.^(node_coords_boundary(:,2));
   
   %%%%%%%%% Compute directivity using spectral relation with the field on
   %%%%%%%%% the obstacle surface

   Directivities= sol_ar.'*phase;

   Dir_star = sol_star.'*phase;
   


% %%%%% incidence parameter to calculate
%  beta_star = 0.8;
%  theta_star = atan(beta_star);
%  [s_star,q_star] = find_plane_wave_parameters(K,beta_star);

 tildeDirectivities = 0*Directivities;

%%%% let's do some interpolation
%%%% and compute coefficients for embedding formula

Coeffs = zeros(N_angles,N_angles);
rhs = zeros(N_angles,1);

for n_cur = 1: N_angles
        s_cur_in = s_ar(n_cur);
        rhs_factor = (s_star + 1/s_star- s_cur_in - 1/s_cur_in);
 %       abs(rhs_factor)
        rhs_Dir = interp1(theta_ar_obs,Directivities(n_cur,:),theta_star);
        rhs(n_cur) = -rhs_factor*rhs_Dir;
    for m_cur = 1:N_angles
        s_cur = s_ar(m_cur);
        cur_factor = (s_cur + 1/s_cur- s_cur_in - 1/s_cur_in);
        cur_Dir = interp1(theta_ar_obs,Directivities(n_cur,:),theta_ar_in(m_cur));
        %%%%%%% coefficients for embedding
        Coeffs(m_cur,n_cur) = cur_factor*cur_Dir;
    end
    tildeDirectivities(n_cur,:) = (s_ar_obs + 1./s_ar_obs - s_ar(n_cur) - 1/s_ar(n_cur)).*Directivities(n_cur,:);
end



factor = (s_ar_obs + 1./s_ar_obs - s_star - 1/s_star);
tildeDir_star = (s_ar_obs + 1./s_ar_obs - s_star - 1/s_star).*Dir_star;

Coeffs_embed = Coeffs\rhs;


%%%%%%%%%% Embedding formula
tildeDir_star_embed = (ones(1,N_angles)*(tildeDirectivities.*repmat(Coeffs_embed,1,N_angles_obs)));

%%%%% Check embedding formula visually

  fig = figure;
   plot(beta_ar_obs,real(tildeDir_star_embed),'*')
   hold on
   plot(beta_ar_obs,real(tildeDir_star),'-')
   plot(beta_ar_obs,imag(tildeDir_star_embed),'*')
   plot(beta_ar_obs,imag(tildeDir_star),'-')
   %title("Directivity",FontSize=18)
   xlabel('$\beta$',FontSize=16,Interpreter='latex')
   ylabel('${\rm Re}[\tilde S(\beta,\beta^{\rm in})]$, ${\rm Im}[\tilde S(\beta,\beta^{\rm in})]$',FontSize=16,Interpreter='latex')
   axis([0 4.5 -inf inf])
   legend('${\rm Re}[\tilde S(\beta,\beta^{\rm in})]$ from embedding',...
       '${\rm Re}[\tilde S(\beta,\beta^{\rm in})]$ from BAE',...
       '${\rm Im}[\tilde S(\beta,\beta^{\rm in})]$ from embedding',...
       '${\rm Im}[\tilde S(\beta,\beta^{\rm in})]$ from BAE', ... 
       FontSize=16,Interpreter='latex')
 %pause
    % exportgraphics(fig, 'directivity.pdf', ...
    % 'ContentType', 'vector')
% vq = interp1(theta_ar_obs,Directivities(1,:),theta_star);
% 
% 
% hold all
% 
% plot(theta_star,abs(vq),'*')



