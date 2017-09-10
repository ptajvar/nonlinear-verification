% load PH
close all
PH = Polyhedron('V',[-1 -1;-1 1;1 -1;1 1]);
pwa = cell(1);
pwa{1}.nx = 2;
pwa{1}.nu = 2;
pwa{1}.grid = {[-1 1],[-1 1],[-.5 .5],[-.5 .5]};
pwa{1}.map = {[-1],[-1],[-.5],[-.5]};
pwa{1}.ulim = [-.5 .5;-.5 .5];
pwa{1}.A = 2*eye(2);
pwa{1}.C = zeros(0,2);
pwa{1}.D = zeros(0,2);
pwa{1}.Bu = eye(2);
pwa{1}.f = zeros(2,1);
pwa{1}.Bw = 0*[-1 1;-1 1];
H = PH.H;
SpecA = H(:,1:end-1);
SpecB = H(:,end);
subplot(1,2,1);
plot(PH);
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
title('(a)')
xlabel('x_1')
ylabel('x_2')

% second system
subplot(1,2,2);
pwa2 = pwa;
pwa2{1}.Bw = .2*[-1 1;-1 1];
plot(PH);
hold on
problem = safe_control_box(SpecA,SpecB,pwa2);
solution = problem.getSol;
for i = 1:length(solution)
    res(i) = solution{i,2};
end
Mu = res(2:3)';
G = res(4:5)';
resPH = Polyhedron('A',[eye(2);-eye(2)],'b',[Mu+G;-(Mu-G)]);
plot(resPH,'color','green');
title('(b)')
xlabel('x_1')
ylabel('x_2')
legend('Specification','Invariant Box')

% 
% 
% [ controlled_sys, controller ] = Add_Controller(pwa2,solution);
% 
% [X,Y] = ginput(1);
% state = [X;Y];
% %% Plotting the evolution of the state
% nIter = 20;
% for i = 1:nIter
%     realState = eval_pwa(controlled_sys,state(:,end));
%     state = [state realState];
%     plot(realState(1,:),realState(2,:),'Color',[ones(1,3)*i/nIter],'Marker','*');
end