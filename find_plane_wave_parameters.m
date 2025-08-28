function [s_in,q_in] = find_plane_wave_parameters(K,beta)

[s1_in,s2_in,s3_in,s4_in] = saddle_points_torus(K,beta);

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
end