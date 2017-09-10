function [next,err,Y] = eval_pwa(pwa,X,U)
%% Init
nx = pwa{1}.nx;
nu = pwa{1}.nu;
if nargin<3 || isempty(U)
    U = zeros(nu,size(X,2));
end

if size(X,1)~=nx || size(U,1)~=nu
    error('Dimension mismatch between PWA system and states/inputs');
end
Xindex = cell(1,nx); % auxilary variable to map state values to partition subscripts
Uindex = cell(1,nu); % auxilary variable to map input values to partition subscripts
next = zeros(size(X));
err = zeros(size(X));
%% Compute the next state for each state/input pair
for it = 1:size(X,2) % each column of X corresponds to a state vector
    %% finding the corresponding partition of the input point
    for ii = 1:nx
        Xindex{ii} = find(pwa{1}.map{ii}<=X(ii,it),1,'last');
        if (isempty(Xindex{ii}))
            Xindex{ii} = 1;
            warning(['Point X(' num2str(ii) ') = ' X(ii,it) ' out of grid bound']);
        end
    end
    for ii = 1:nu
        Uindex{ii} = find(pwa{1}.map{nx+ii}<=U(ii,it),1,'last');
        if (isempty(Uindex{ii}))
            Uindex{ii} = 1;
            warning(['Point X(' num2str(ii) ') = ' U(ii,it) ' out of grid bound']);
        end
    end
    %% Compute the affine function
    vindex = [Xindex, Uindex];
    err = pwa{vindex{:}}.Bw;
    nextTemp = pwa{vindex{:}}.f+pwa{vindex{:}}.A*X(:,it)+...
        pwa{vindex{:}}.Bu*U(:,it);
    next(:,it) = nextTemp;    
    Y(:,it) = pwa{vindex{:}}.C*nextTemp;
end