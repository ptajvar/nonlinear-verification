f = delay_logistics_gen;
ff = {@(X,Y) f(X,Y), @(X,Y) X};
C = [1 0];
D = zeros(1,0 + 2 + 1); % p x (m + n + 1)
xmin = 1;
xmax = 2;
step = 0.01;
meshSize = [3 3];
xv{1} = linspace(xmin,xmax,meshSize(1));
xv{2} = linspace(xmin,xmax,meshSize(2));
uv = {};
pwa_container_p = pwa_estimator( ff, xv, uv, C, D, 'p' );
pwa_container_l = pwa_estimator( ff, xv, uv, C, D, 'l' );
%% drawing main nonlinear function
figure,
x = xmin:step:xmax;
[X_N_1,X_N_2] = meshgrid(x);
mesh(X_N_1,X_N_2,f(X_N_1,X_N_2));
title('nonlinear function');
xlabel('X^{n-1}');
ylabel('X^{n-2}');
zlabel('X^{n}');
%% drawing projection pwa function with error bounds
figure,
subplot(1,2,1);
[X,E,Y] = eval_pwa(pwa_container_p,[X_N_1(:);X_N_2(:)],[]);
mesh(X_N_1,X_N_2,Y+E);
hold on
mesh(X_N_1,X_N_2,Y-E);
title('projection pwa function')
xlabel('X^{n-1}');
ylabel('X^{n-2}');
zlabel('X^{n}');
%% drawing lipschitz pwa function with error bounds
subplot(1,2,2);
[X,E] = eval_pwa(pwa_container_l,X_N_1,X_N_2);
mesh(X_N_1,X_N_2,X+E);
hold on
mesh(X_N_1,X_N_2,X-E);
title('lipschitz pwa function')
xlabel('X^{n-1}');
ylabel('X^{n-2}');
zlabel('X^{n}');