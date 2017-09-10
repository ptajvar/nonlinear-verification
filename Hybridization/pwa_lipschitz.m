function pwa_container = pwa_lipschitz(fhandle,xv,uv)

% Inputs,
% - fhandle: handle to the nonlinear function to be estimated
% - xv: cell array (1,nx) where each cell contains a vector defining the
% partiotion intervals on that state
% - uv: cell array (1,nu) where each cell contains a vector defining the
% partiotion intervals on that input
%
% Outputs,
% - estimation: function handle to the pwa estimation
% - pwa_container: cell array of the pwa estimation (the format suitable
% for passing to 'writePWA.m')

% todo: insert pre condition verification size(range) = (length(grid) X 2)
% todo: insert pre condition min(gridSize)>1

fv = [xv uv];
nx = length(xv);
nu = length(uv);
dim = length(fv);
range = zeros(dim,2);
gridSize = zeros(1,dim);
for i = 1:dim
    gridSize(i) = length(fv{i});
end
%% initializing container
pwa_container = cell((gridSize-1)); % cell container to be passed to 'writePWA.m'
pwa_container{1}.nu = length(uv); % number of inputs
pwa_container{1}.nx = length(xv); % number of states
pwa_container{1}.ny = 1; % number of outputs
pwa_container{1}.xmin = zeros(length(xv));
pwa_container{1}.xmax = zeros(length(xv));
pwa_container{1}.umin = zeros(length(uv));
pwa_container{1}.umax = zeros(length(uv));
for i = 1:dim
    range(i,1) = fv{i}(1);
    range(i,2) = fv{i}(end);
    gridSize(i) = length(fv{i});
    if i<=length(xv)
        pwa_container{1}.xmin(i) = range(i,1);
        pwa_container{1}.xmax(i) = range(i,2);
    else
        pwa_container{1}.umin(i-length(xv)) = range(i,1);
        pwa_container{1}.umax(i-length(xv)) = range(i,2);
    end
end
X = inf*ones(dim,max(gridSize)); % marks the start of each partition (region/piece) so starts from -inf to mark out of bound region where the closest region rules shall be applied
for i = 1:dim
    X(i,1:gridSize(i)) = fv{i}(:);
    pwa_container{1}.map{i} = fv{i}(1:end-1); % a map that is used to refer a point to its corresponding partition (region)
    pwa_container{1}.grid{i} = fv{i};
end
subscript = cell(1,dim); % auxilary cell to allow iterating on dim-D matrices

%% calculation of output value for each piece of grid (middle point of the hyper cube)
% aka creating pwa estimation
sample = struct('A', zeros(1,dim), 'const',0);
S = gridSize-1; % size of mesh of regions
result = repmat(sample,S);
leftVal = zeros(dim,1); % vector of lesser values in a region
rightVal = zeros(dim,1);% vector of higher values in a region
for i = 1:numel(result)
    [subscript{:}] = ind2sub(S,i);
    for j = 1:dim
        leftVal(j) = X(j,subscript{j});
        rightVal(j) = X(j,subscript{j}+1);
    end
    midVal = (leftVal+rightVal)/2;
    result(i).const = fhandle(midVal(1:nx),midVal(nx+1:end));
    for j = 1:dim
        lower = midVal;
        higher = midVal;
        lower(j) = leftVal(j);
        higher(j) = rightVal(j);
        L(j) = higher(j)-lower(j);
        lowVal = fhandle(lower(1:nx),lower(nx+1:end));
        highVal = fhandle(higher(1:nx),higher(nx+1:end));
        result(i).A(j) = ((highVal-lowVal)/L(j));
    end
    result(i).const = result(i).const - result(i).A*midVal;
    w = [0 0]; % maximum estimation error value (to be calculated through estimation)
    nPoints = 100;
    r = rand(dim,nPoints);
    d1 = ceil(rand(1,nPoints)*dim);
    d2 = ceil(rand(1,nPoints)*dim);
    diffRat = 0.01; % ratio of distance between two nearby points for evaluation of gradient to the length of the mode
    MH = 0;
    mH = 0;
    for j = 1:nPoints
        point = (1-r(:,j)).*leftVal+r(:,j).*rightVal;
        expected = result(i).const+result(i).A*point;
        actual = fhandle(point(1:nx),point(nx+1:end));
        error = actual-expected;
        w(1) = min(w(1),error);
        w(2) = max(w(2),error);
%         diff1 = zeros(dim,1);
%         diff1(d1(j)) = L(d1(j))*diffRat;
%         diff2 = zeros(dim,1);
%         diff2(d2(j)) = L(d2(j))*diffRat;
%         point2 = point+diff1;
%         point3 = point+diff2;
%         point4 = point3+diff1;
%         tan1 = (fhandle(point2(1:nx),point2(nx+1:end))-fhandle(point(1:nx),point(nx+1:end)))...
%             /max(diff1);
%         tan2 = (fhandle(point4(1:nx),point4(nx+1:end))-fhandle(point3(1:nx),point3(nx+1:end)))...
%             /max(diff1);
%         lap = (tan2-tan1)/max(diff2);
%         MH = max(MH,lap);
%         mH = min(mH,lap);
    end
%     w = [0 0];
%     for j = 1:dim
%         for k = 1:dim
%             w(2) = w(2) + 1/8*MH*L(j)*L(k);
%             w(1) = w(1) + 1/8*mH*L(j)*L(k);
%         end
%     end  
    pwa_container{i}.A = result(i).A(1:nx);
    pwa_container{i}.Bu = result(i).A(nx+1:end);
    pwa_container{i}.f = result(i).const;
    pwa_container{i}.Bw = w*1.2;
    pwa_container{i}.lim = [leftVal rightVal];
end