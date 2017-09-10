% load PH
close all
xmin = 1;
xmax = 2;
f = @(x,u) x.^3+u;
D = zeros(2);
C = eye(2);
PH = Polyhedron('V',[xmin xmin;xmin xmax;xmax xmin;xmax xmax]);
ulim = 1*[-1 1;-1 1];
nmod = [2 2];
Xgrid = {linspace(xmin,xmax,nmod(1)+1),linspace(xmin,xmax,nmod(2)+1)};
Ugrid = {ulim(1,:),ulim(2,:)};
pwa = pwa_estimator( f, Xgrid, Ugrid, C, D, 'l' );
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

% legend('Specification','Invariant Box')

% 
% [ controlled_sys, controller ] = Add_Controller(pwa,solution);
% 
% [X,Y] = ginput(1);
% state = [X;Y];
% %% Plotting the evolution of the state
% nIter = 20;
% for i = 1:nIter
%     input(1:2,i) = eval_pwa(controller,[0;0],state(:,end));
%     state = [state f(state(:,end),input(1:2,i))];    
%     realState = state(:,end-1);
%     plot(realState(1,:),realState(2,:),'Color',[ones(1,3)*i/nIter],'Marker','*');
% end

nmod = [10 10];
Xgrid = {linspace(xmin,xmax,nmod(1)+1),linspace(xmin,xmax,nmod(2)+1)};
Ugrid = {ulim(1,:),ulim(2,:)};
pwa = pwa_estimator( f, Xgrid, Ugrid, C, D, 'l' );
pwa{1}.ulim = ulim;
H = PH.H;
SpecA = H(:,1:end-1);
SpecB = H(:,end);
subplot(1,2,2);
plot(PH);
title('(b)')
hold on
problem = safe_control_box(SpecA,SpecB,pwa);
solution = problem.getSol;
clear res
for i = 1:5
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