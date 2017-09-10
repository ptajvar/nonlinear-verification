% PH = polyhPoints;
% save PH PH
load PH
figure,
H = PH.H;
SpecA = H(:,1:end-1);
SpecB = H(:,end);
grid{1} = [-10 10];
grid{2} = [-10 10];
pwa{1}.grid = grid;
pwa{1}.A = eye(2);
pwa{1}.f = zeros(2,1);
pwa{1}.Bu = zeros(2,0);
pwa{1}.C = zeros(0,2);
pwa{1}.D = zeros(0,0);
pwa{1}.nx = 2;
pwa{1}.nu = 0;
pwa{1}.ulim = zeros(0,2);
pwa{1}.Bw = zeros(2,2);


subplot(1,2,1)
plot(PH);
title('Before Rotation');
hold on
problem1 = safe_control_box(SpecA,SpecB,pwa);
solution = problem1.getSol;
for i = 1:length(solution)
    res(i) = solution{i,2};
end
Mu = res(2:3)';
G = res(4:5)';
resPH = Polyhedron('A',[eye(2);-eye(2)],'b',[Mu+G;-(Mu-G)]);
plot(resPH,'color','green');
legend('Specification','Inner Box')

[T,dc] = SpacePCA(SpecA,SpecB);
SpecArot = SpecA*T^-1;
SpecBrot = SpecB;  %  - SpecA*dc'
pwa2 = pwa;
pwa2{1}.A = T^-1*pwa{1}.A*T;
problem2 = safe_control_box(SpecArot,SpecBrot,pwa2);
solution = problem2.getSol;
for i = 1:length(solution)
    res2(i) = solution{i,2};
end
Mu2 = res2(2:3)';
G2 = res2(4:5)';
res2PH = Polyhedron('A',[eye(2);-eye(2)]*T,'b',[Mu2+G2;-(Mu2-G2)]); % +[eye(2);-eye(2)]*T*dc'
subplot(1,2,2)
plot(PH)
title('After Rotation');
hold on
plot(res2PH,'color','green');
legend('Specification','Inner Box')

print(gcf,'Rotation.png','-dpng','-r350')