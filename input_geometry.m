

%%%% square

N_p = 20;

first_side = [-N_p:N_p; -N_p*ones(1,2*N_p+1)];
second_side = [N_p*ones(1,2*N_p+1);-N_p:N_p];
third_side = [N_p:-1:-N_p;N_p*ones(1,2*N_p+1)];
fourth_side = [-N_p*ones(1,2*N_p+1);N_p:-1:-N_p];

square  = [first_side,second_side,third_side,fourth_side].';

% figure;
% plot(square_boundary(:,1),square_boundary(:,2),'*')

%%%%%% strip

strip = first_side.';

%strip_boundary = second_side.';

%%%%%% P-shaped figure
first_side_p = [-5:5; -5*ones(1,2*5+1)];
second_side_p = [5*ones(1,2*5+1);-5:5];
third_side_p = first_side_p+[10;10];
fourth_side_p = [16*ones(1,22);5:-1:-16];
fifth_side_p = [16:-1:-5;-16*ones(1,22)];
sixth_side_p = second_side_p + [-10;-10];

P_figure = [first_side_p,second_side_p,third_side_p,fourth_side_p,...
    fifth_side_p,sixth_side_p].';
 % figure;
 % plot(P_figure(:,1),P_figure(:,2),'*')

%%%%% G-shaped strip

g_strip = unique([first_side,second_side].','rows');

  % figure;
  % plot(g_strip(:,1),g_strip(:,2),'*')