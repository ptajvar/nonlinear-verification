% load PH
close all
PH = Polyhedron('V',[0 0;0 .8;1 0;1 1]);
pwa = cell(2,2);
pwa{1}.grid = {[0 .5 1],[0 .5 1]};
pwa{1}.ulim = 1*[-1 1;-1 1];
pwa{1}.Bu = eye(2);
pwa{2}.Bu = eye(2);
pwa{3}.Bu = eye(2);
pwa{4}.Bu = eye(2);
Wlim = .2*[-1 1;-1 1];
pwa{1}.Wlim = Wlim;
pwa{2}.Wlim = Wlim;
pwa{3}.Wlim = Wlim;
pwa{4}.Wlim = Wlim;
pwa2 = pwa;
for i = 1:numel(pwa)
    pwa{i}.f = zeros(2,1);
    pwa2{i}.f = zeros(2,1);
    pwa{i}.Bw = zeros(2,2);
    pwa2{i}.Bw = zeros(2,2);
    pwa2{i}.A = 1*eye(2);
end
pwa = pwa2;
pwa{2}.A = [3 0;0 3];
pwa{1}.A = [30 0;0 30];
pwa{4}.A = [30 0;0 30];
pwa{3}.A = [30 0;0 30];
H = PH.H;
SpecA = H(:,1:end-1);
SpecB = H(:,end);

subplot(1,2,1)
plot(PH);
title('Without Dybamics');
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
legend('Specification','Invariant Box')

subplot(1,2,2)
plot(PH)
title('With Dynamics');
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
legend('Specification','Invariant Box')


print(gcf,'Dynamics.png','-dpng','-r500');