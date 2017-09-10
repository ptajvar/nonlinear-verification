function prob = safe_control_box(SpecA, SpecB, pwa)

% Require SpecA \in R^{ns * nd} and SpecB \in R^{ns} where 'ns' is number
% of inequality specifications and nd is the number of space dimensions.
% grid is a set of adjacent interval in each dimension of the state space:
% I_i = {[x_i_min,x_i_1),[x_i_1,x_i_2),...[x_i_nd-1,x_i_max),[x_i_max,x_i_max]}
% Output MILP is a structure that holds a mixed-integer linear program that
% finds the largest box inside the Specification and marks the sections of
% the grid that have nonempty intersection with the Box.

P = MILPC; % P is an object that holds MILP variables and constraints

grid = pwa{1}.grid;
SpecSize = size(SpecA);
nc = SpecSize(1); % number of specification inequalities
nx = SpecSize(2); % number of space dimensions
nu = size(pwa{1}.Bu,2);
PWASize = size(pwa);
grid(nx+1:end) = [];
PWASize(nx+1:end) = [];
if nx ~= length(grid)
    error('Grid definition does not match space dimension');
else
    nbin = zeros(1,nx); % number of binary variables corresponding to
    % intersection of an interval with the box
    loBound = zeros(nx,1);
    hiBound = zeros(nx,1);
    for i = 1:nx
        nbin(i) = length(grid{i})-1;
        loBound(i) = grid{i}(1);
        hiBound(i) = grid{i}(end);
    end
end
csumBin = cumsum(nbin); % to find the interval to which the variable belongs
prodBin = prod(nbin);
%% Defining the Box
P = P.Add_Variable('SmallestGenerator','double',1,[0 min((hiBound-loBound)/2)]);
P = P.Add_Variable('Mu','double',nx,[loBound,hiBound]); % Mu is the center of the box
P = P.Add_Variable('G','double',nx,[zeros(nx,1), (hiBound-loBound)/2]); % g is the generator length
%% Cost function f:
P = P.Set_Cost({'SmallestGenerator',{'G',1:nx}},[-1 -.01*ones(1,nx)]);
%% Building Constraint matrices A and b
% Part 1: SmallestGenerator<=g(i)
for i = 1:nx
    P = P.Add_Constraint({'SmallestGenerator', {'G',i}},[1 -1],0);
end
% Part 2a: Box is inside Specification: SpecA*x<SpecB \forall x \in Box
% for each constraint we only need to verify the corner that maximises
% SpecA(i,:)*x
for i = 1:nc
    P = P.Add_Constraint({{'Mu',1:nx},{'G',1:nx}},[SpecA(i,:) abs(SpecA(i,:))],SpecB(i));
end
% Part 2b: Box is inside the modes
for i = 1:nx
    P = P.Add_Constraint({{'Mu',i},{'G',i}},[1 1],hiBound(i));
    P = P.Add_Constraint({{'Mu',i},{'G',i}},[-1 1],-loBound(i));
end
% Part 3: Finding the intervals that have nonempty intersection with the
% box (Big M method)
M = 1e3;
for dim = 1:nx
    name_dlo = ['delta_lo ' num2str(dim)]; % name of variables in this dimension
    name_dhi = ['delta_hi ' num2str(dim)];
    name_lo = ['LOWER ' num2str(dim)];
    name_hi = ['HIGHER ' num2str(dim)];
    P = P.Add_Variable(name_dlo,'bin',nbin(dim)+1); % +1 because number of end points in each dimension is (number_of_intervals+1)
    P = P.Add_Variable(name_dhi,'bin',nbin(dim)+1); % +1 because number of end points in each dimension is (number_of_intervals+1)
    P = P.Add_Variable(name_lo,'double',nbin(dim)); % referred to as \alpha in invariant box chapter
    P = P.Add_Variable(name_hi,'double',nbin(dim)); % referred to as \beta in invariant box chapter
    for j = 1:nbin(dim)+1
        % for the interval to intersect with, two conditions are necessary:
        % 1- The beginning of the box should be before end of the interval (delta_lo)
        % 2- The beginning of the interval should be before end of the box (delta_hi)
        P = P.Add_Constraint({{name_dlo,j},{'Mu',dim},{'G',dim}},...
            [M 1 -1],grid{dim}(j)+M);
        P = P.Add_Constraint({{name_dlo,j},{'Mu',dim},{'G',dim}},...
            [-M -1 1],-grid{dim}(j));
        P = P.Add_Constraint({{name_dhi,j},{'Mu',dim},{'G',dim}},...
            [M -1 -1],-grid{dim}(j)+M);
        P = P.Add_Constraint({{name_dhi,j},{'Mu',dim},{'G',dim}},...
            [-M 1 1],grid{dim}(j));
        % Introducing lower/higher bounds of the box
        if j>nbin(dim) 
            continue;
        end
        P = P.Add_Constraint({{name_dlo,j},{name_lo,j}},...
            [M 1],M+grid{dim}(j));
        P = P.Add_Constraint({{name_dlo,j},{name_lo,j}},...
            [M -1],M-grid{dim}(j));
        P = P.Add_Constraint({{name_dlo,j},{name_lo,j},{'Mu',dim},{'G',dim}},...
            [-M 1 -1 1],0);
        P = P.Add_Constraint({{name_dlo,j},{name_lo,j},{'Mu',dim},{'G',dim}},...
            [-M -1 1 -1],0);
        
        P = P.Add_Constraint({{name_dhi,j+1},{name_hi,j}},...
            [M 1],M+grid{dim}(j+1));
        P = P.Add_Constraint({{name_dhi,j+1},{name_hi,j}},...
            [M -1],M-grid{dim}(j+1));
        P = P.Add_Constraint({{name_dhi,j+1},{name_hi,j},{'Mu',dim},{'G',dim}},...
            [-M 1 -1 -1],0);
        P = P.Add_Constraint({{name_dhi,j+1},{name_hi,j},{'Mu',dim},{'G',dim}},...
            [-M -1 1 1],0);
    end
end

%% Dynamics Verification
ulim = pwa{1}.ulim;
for i = 1:prodBin
    name = ['fc_' num2str(i)];
    P = P.Add_Variable(name,'double',nu);
    name = ['absfc_' num2str(i)]; % absolute valu of fc
    P = P.Add_Variable(name,'double',nu);
    for j = 1:nu
        name = ['Ac_' num2str(i) '_' num2str(j)];
        P = P.Add_Variable(name,'double',nx);
        % ulim constraint (u \in U)
        P = P.Add_Constraint({{['absfc_' num2str(i)],j},...
            {['fc_' num2str(i)],j}},[-1 1],eps);
        P = P.Add_Constraint({{['absfc_' num2str(i)],j},...
            {['fc_' num2str(i)],j}},[-1 -1],eps);
        name = ['absAc_' num2str(i) '_' num2str(j)]; % absolute valu of Ac
        P = P.Add_Variable(name,'double',nx);
        for it = 1:nx
            P = P.Add_Constraint({{['absAc_' num2str(i) '_' num2str(j)],it},...
                {['Ac_' num2str(i) '_' num2str(j)],it}},[-1 1],0);
            P = P.Add_Constraint({{['absAc_' num2str(i) '_' num2str(j)],it},...
                {['Ac_' num2str(i) '_' num2str(j)],it}},[-1 -1],0);
        end
        P = P.Add_Constraint({{['absAc_' num2str(i) '_' num2str(j)],1:nx},...
            {['fc_' num2str(i)],j}},[ones(1,nx)/2 1],ulim(j,2));
        P = P.Add_Constraint({{['absAc_' num2str(i) '_' num2str(j)],1:nx},...
            {['fc_' num2str(i)],j}},[ones(1,nx)/2 -1],-ulim(j,1));
        %end - ulim constraint
    end
    % involved deltas' indices will be stored in k
    [I{1:nx}] = ind2sub(PWASize,i);
    % dynamical constraints
    A = pwa{i}.A;
    f = pwa{i}.f;
    B = pwa{i}.Bu;
    Wlim = [-pwa{i}.Bw pwa{i}.Bw];
    
    % making abs((A diag(p)-B*Ac))
    for j = 1:nx % row of A
        name = ['absA_' num2str(i) '_' num2str(j)];
        P = P.Add_Variable(name,'double',nx);
        for jj = 1:nx % column of A
            name_lo = ['LOWER ' num2str(jj)];
            name_hi = ['HIGHER ' num2str(jj)];
            clear vars multipliers
            vars{1} = {name,jj};
            multipliers(1) = -1;
            vars{2} = {name_lo,I{jj}};
            multipliers(2) = -A(j,jj)/2;
            vars{3} = {name_hi,I{jj}};
            multipliers(3) = A(j,jj)/2;
            for iu = 1:nu
                vars{end+1} = {['Ac_' num2str(i) '_' num2str(iu)],jj};
                multipliers(end+1) = B(j,iu);
            end
            P = P.Add_Constraint(vars,multipliers,0);
            multipliers(2:end) = -multipliers(2:end);
            P = P.Add_Constraint(vars,multipliers,0);
        end
    end
    % END - making abs((A-B*Ac))*(hi-lo)
    % Verifying box invariance
    for j = 1:nx
        A_dim = A(j,:);
        clear vars multipliers
        for jj = 1:nx
            name_lo = ['LOWER ' num2str(jj)];
            name_hi = ['HIGHER ' num2str(jj)];
            name_dlo = ['delta_lo ' num2str(jj)];
            name_dhi = ['delta_hi ' num2str(jj)];
            vars{jj} = {name_lo,I{jj}};
            multipliers(jj) = A_dim(jj)/2;
            vars{jj+nx} = {name_hi,I{jj}};
            multipliers(jj+nx) = A_dim(jj)/2;
            vars{jj+2*nx} = {name_dlo,I{jj}+1};
            multipliers(jj+2*nx) = M;
            vars{jj+3*nx} = {name_dhi,I{jj}};
            multipliers(jj+3*nx) = M;
        end
        vars{end+1} = {['absA_' num2str(i) '_' num2str(j)],1:nx};
        multipliers = [multipliers ones(1,nx)];
        vars{end+1} = {['fc_' num2str(i)],1:nu};
        multipliers = [multipliers B(j,:)];
        vars{end+1} = {'Mu',j};
        multipliers = [multipliers -1];
        vars{end+1} = {'G',j};
        multipliers = [multipliers -1];
        P = P.Add_Constraint(vars,multipliers,M*2*nx+eps-Wlim(j,2)-f(j));
        
        clear vars multipliers
        for jj = 1:nx
            name_lo = ['LOWER ' num2str(jj)];
            name_hi = ['HIGHER ' num2str(jj)];
            name_dlo = ['delta_lo ' num2str(jj)];
            name_dhi = ['delta_hi ' num2str(jj)];
            vars{jj} = {name_lo,I{jj}};
            multipliers(jj) = -A_dim(jj)/2;
            vars{jj+nx} = {name_hi,I{jj}};
            multipliers(jj+nx) = -A_dim(jj)/2;
            vars{jj+2*nx} = {name_dlo,I{jj}+1};
            multipliers(jj+2*nx) = M;
            vars{jj+3*nx} = {name_dhi,I{jj}};
            multipliers(jj+3*nx) = M;
        end
        vars{end+1} = {['absA_' num2str(i) '_' num2str(j)],1:nx};
        multipliers = [multipliers ones(1,nx)];
        vars{end+1} = {['fc_' num2str(i)],1:nu};
        multipliers = [multipliers -B(j,:)];
        vars{end+1} = {'Mu',j};
        multipliers = [multipliers 1];
        vars{end+1} = {'G',j};
        multipliers = [multipliers -1];
        P = P.Add_Constraint(vars,multipliers,M*2*nx+eps+Wlim(j,1)+f(j));
    end
    % END-Verifying box invariance
end

% MILP = P.getMILP;
prob = P;