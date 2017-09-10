% load PH
clear
close all
xmin = [-6 -10];
xmax = [8 10];
f = rakovic1;
D = zeros(2);
C = eye(2);
ulim = [-1 1];
nmod = [2 2];

Xgrid = {[-6 1 8],[-10 10]};
Ugrid = {[ulim(1),ulim(2)]};
pwa = pwa_estimator( f, Xgrid, Ugrid, C, D, 'l' );

pwa{1}.ulim = ulim;
for i = 1:2
    pwa{i}.Bw = [-.1 .1;-.1 .1];
end
SpecA = [-1 1;-3 -1;.2 1;-1 0;1 0];
SpecB = [15;25;9;6;8];
H = Polyhedron('A',SpecA,'b',SpecB);
plot(H);
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
axis([-6 8 -10 10]);
legend('Specification','Invariant Box')

[ controlled_sys, controller ] = Add_Controller(pwa,solution);
[X,Y] = ginput(1);
state = [X;Y];
%% Plotting the evolution of the state
nIter = 15;
for i = 1:nIter
    state = [state eval_pwa(controlled_sys,state(:,end))];
    realState = state(:,end-1);
    plot(realState(1,:),realState(2,:),'Color',[ones(1,3)*i/nIter],'Marker','*');
end