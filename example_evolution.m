clear
f = delay_logistics_gen;
C = [1 0];
D = zeros(1,0 + 2 + 1); % p x (m + n + 1)
xmin = 1.4;
xmax = 2;
step = 0.01;
meshSize = [6 6];
xv{1} = linspace(xmin,xmax,meshSize(1));
xv{2} = linspace(xmin,xmax,meshSize(2));
uv = {};
pwa_container = pwa_estimator( f, xv, uv, C, D, 'l' );
%% drawing main nonlinear function
T = 100;
x = zeros(1,T);
x(1) = 1.7;
x(2) = 1.9;
for i = 3:T
    x(i) = [1,0]*f([x(i-1);x(i-2)]);
end
plot (1:T,x,'LineWidth',1)
xlabel('time k');
ylabel('x(k)');
hold on
%% drawing pwa function with error bounds
xmin = nan(1,T);
xmax = nan(1,T);
xmin(1) = 1.7;
xmax(2) = 1.9;
for it = 1:1000
    for i = 3:T
        r = rand*2-1;
        [X,E] = eval_pwa(pwa_container,[x(i-1);x(i-2)]);
        x(i) = X(1) + r*E(1);
    end
    xmin = min(x,xmin);
    xmax = max(x,xmax);
end
plot (1:T,xmin,'LineWidth',1)
plot (1:T,xmax,'LineWidth',1)
legend('nonlinear function','pwa minimum bound','pwa maximum bound')