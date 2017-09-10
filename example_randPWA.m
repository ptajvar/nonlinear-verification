% load PH
close all
xmin = -1;
xmax = 1;
f = @(x,u) x.^3+u;
D = zeros(2);
C = eye(2);
PH = Polyhedron('V',[xmin xmin;xmin .7*xmax;xmax 0.5*xmin;xmax xmax]);
ulim = 1*[-1 1;-1 1];
nmod = [2 2];
Xgrid = {linspace(xmin,xmax,nmod(1)+1),linspace(xmin,xmax,nmod(2)+1)};
Ugrid = {ulim(1,:),ulim(2,:)};
pwa = pwa_estimator( f, Xgrid, Ugrid, C, D, 'l' );

% define A matrix
pwa{1}.A = 3*eye(2);
pwa{2}.A = eye(2);
pwa{3}.A = eye(2);
pwa{4}.A = eye(2);
% end = define A matrix
pwa{1}.ulim = ulim;
H = PH.H;
SpecA = H(:,1:end-1);
SpecB = H(:,end);
subplot(1,2,1);
plot(PH);
title('(a)')
hold on
problem = safe_control_box(SpecA,SpecB,pwa);
solution = problem.getSol;
for i = 1:length(solution)
    res(i) = solution{i,2};
end
Mu = res(2:3)';
G = res(4:5)';
resPH = Polyhedron('A',[eye(2);-eye(2)],'b',[Mu+G;-(Mu-G)]);
plot(resPH,'color','green');
xlabel('x_1')
ylabel('x_2')

for i = 1:length(Xgrid{1})
    plot([Xgrid{1}(i) Xgrid{1}(i)],[Xgrid{2}(1) Xgrid{2}(end)],'k');
end

for i = 1:length(Xgrid{2})
    plot([Xgrid{1}(1) Xgrid{1}(end)],[Xgrid{2}(i) Xgrid{2}(i)],'k');
end


Xgrid = {linspace(xmin,xmax,nmod(1)+1),linspace(xmin,xmax,nmod(2)+1)};
Ugrid = {ulim(1,:),ulim(2,:)};

pwa{1}.ulim = ulim;
% define A matrix
pwa{1}.A = 30*eye(2);
pwa{2}.A = eye(2);
pwa{3}.A = eye(2);
pwa{4}.A = eye(2);
% end = define A matrix
H = PH.H;
SpecA = H(:,1:end-1);
SpecB = H(:,end);
subplot(1,2,2);
plot(PH);
title('(b)')
hold on
problem = safe_control_box(SpecA,SpecB,pwa);
solution = problem.getSol;
for i = 1:length(solution)
    res(i) = solution{i,2};
end
Mu = res(2:3)';
G = res(4:5)';
resPH = Polyhedron('A',[eye(2);-eye(2)],'b',[Mu+G;-(Mu-G)]);
plot(resPH,'color','green');
xlabel('x_1')
ylabel('x_2')

for i = 1:length(Xgrid{1})
    plot([Xgrid{1}(i) Xgrid{1}(i)],[Xgrid{2}(1) Xgrid{2}(end)],'k');
end

for i = 1:length(Xgrid{2})
    plot([Xgrid{1}(1) Xgrid{1}(end)],[Xgrid{2}(i) Xgrid{2}(i)],'k');
end

legend('Specification','Invariant Box')