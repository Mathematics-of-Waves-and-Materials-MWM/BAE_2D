function Gmn = Greens_asympt(K,m,n) 
 
[s1,s2,s3,s4] = saddle_points_torus(K,m/n);

s_ar = [s1,s2,s3,s4];

s = s_ar(imag(s_ar)>0&abs(s_ar)<1&abs(s_ar)>0.8);



q =  -(K^2 - 4 + s + 1/s)/2 + sqrt((K^2-4+s + 1/s)^2 - 4)/2;

if abs(q)>1
    q =  1/q;
end


dq = -q*(s - 1/s)/s/(q - 1/q);

d2q = -2*s^(-3)/(1 - q^(-2)) + 2*(1-s^(-2))*q^(-3)*dq/(1 - q^(-2))^2;

d2phi = -m/s^2 - n/q^2*dq + n/q*d2q;


Gmn = s^m*q^n*sqrt(1/d2phi/2/pi)*1/s/(q-1/q);
end