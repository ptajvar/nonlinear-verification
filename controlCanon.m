function [T,At,Bt] = controlCanon(A,B)

C = ctrb(A,B);
AA = (C^-1*A*C)';
Cd = zeros(length(A),1);
Cd(end) = 1;
for i = 1:length(A)-1
    Cd = [Cd AA*Cd(:,end)];
end
T = Cd*C^-1;
At = T*A*T^-1;
Bt = T*B;