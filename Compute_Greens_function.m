%clear all

%%%%% Greens function is evaluated by direct evaluation of the contour
%%%%% integral. The accuracy can be improved by contour deformation in the
%%%%% steepest descent contour
%%%%% Will fail for real K as branch points will hit the contour

m = 101;
n = 101;
dtheta = 0.001;


theta_ar = 0:dtheta:2*pi;

s_ar = exp(1i*theta_ar);
q_ar = -(K^2 - 4 + s_ar + 1./s_ar)/2 + sqrt((K^2-4+s_ar + 1./s_ar).^2 - 4)/2;

mask = abs(q_ar)>1;

q_ar(mask) = -(K^2 - 4 + s_ar(mask) + 1./s_ar(mask))/2 - sqrt((K^2-4+s_ar(mask) + 1./s_ar(mask)).^2 - 4)/2;

ds_ar = s_ar(2:end) - s_ar(1:end-1);

q_ar = average(q_ar);
s_ar = average(s_ar);


 Greens_func = zeros(m);
 Greens_as = zeros(m);

for n_cur = 1 : n
    for m_cur = 1:m
        Greens_func(m_cur,n_cur) = 1/2i/pi*(s_ar.^(m_cur-2).*q_ar.^(n_cur-1)./(q_ar -1./q_ar))*ds_ar.';
    end
end

% save('Green_function.mat','Greens_func','K')

% 
% figure;
% 
% surf(imag(Greens_func))

%% below some simple checks of accuracy of computations



% m_cur = 3;
% n_cur = 4;
% 
% Greens_asympt(K,m_cur,n_cur)
% 
% Greens_func(m_cur+1,n_cur+1)
% 
% (K^2-4)*Greens_func(m_cur,n_cur) + Greens_func(m_cur+1,n_cur)+Greens_func(m_cur-1,n_cur)...
%     +Greens_func(m_cur,n_cur-1) + Greens_func(m_cur,n_cur+1)
% 
% (K^2-4)*Greens_func(1,1) + 2*Greens_func(2,1) + 2*Greens_func(1,2)

% 
% figure;
% plot(s_ar)
% axis equal