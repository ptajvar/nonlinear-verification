clear
rng(1)
figure
% A = 5*rand(2);
% B = rand(2,1);
A = [1 1;0 1];
B = [0;1];
[TT,AA,BB] = controlCanon(A,B);
% TT = [1 0;0 sum(abs(AA(end,:)))]*TT*[1 0;0 sum(abs(AA(end,:)))]^-1;
% TT = eye(2);
f = twoD_lin_gen;
nx = 2;
nu = 1  ;
C = [1 0];
D = zeros(1);
xmin = [-1 -1];
xmax = [1 1];
umin = -1;
umax = 1;
xv = cell(1,nx);
uv = cell(1,nu);
XmeshSize = [2 2];
UmeshSize = [2 2];
Theta = 0;
T = TT^-1;
Tn = T^-1;
xv = cell(1,nx);
for i = 1:nx
    xmin_temp = sum(Tn(i,:).*(Tn(i,:)>0)*xmin(i) + Tn(i,:).*(Tn(i,:)<0)*xmax(i));
    xmax_temp = sum(Tn(i,:).*(Tn(i,:)<0)*xmin(i) + Tn(i,:).*(Tn(i,:)>0)*xmax(i));
    xv{i} = linspace(xmin_temp,xmax_temp,XmeshSize(i));
end
for i = 1:nu
    uv{i} = linspace(umin,umax,UmeshSize(i));
end
ff = @(xT,u) T^(-1)*f(T*xT,u);
pwa_container2 = pwa_estimator( ff, xv, uv, C*T, D, 'l' );
%% drawing main nonlinear function
PH = Polyhedron('V',[-1 -1;-1 1;1 -1;1 1]);
H = PH.H;
SpecA = H(:,1:end-1);
SpecB = H(:,end);
pwa_container2{1}.ulim = [1;1]*[umin,umax];
problem = safe_control_box(SpecA*T,SpecB,pwa_container2);
plot(PH);
xlabel('Position','fontsize',12);
ylabel('Velocity','fontsize',12);
title('');
hold on
for i = 1:length(problem)
    res(i) = problem{i,2};
end
Mu = res(2:3)';
G = res(4:5)';
resPH = Polyhedron('A',[eye(2);-eye(2)]*Tn,'b',[Mu+G;-(Mu-G)]);
plot(resPH,'color','green');
if it == 2
    legend('Specification','Invariant Box')
end
print(gcf,'CanonicalController.png','-dpng','-r350')