mr = 95;
nr = 0;

Dapr =[];
beta_apr = [];
for nr = 0: 95

idx = node_coords_outer(:,1)==nr&node_coords_outer(:,2)==mr;
beta_apr = [beta_apr,nr/mr];
Dapr = [Dapr,scat_field(idx)/Greens_func(mr+1,nr+1)];
end

figure; 
plot(beta_apr,abs(Dapr))
hold all
plot(beta_ar_obs, abs(Dir_star))


