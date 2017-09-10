function [ post_pwa ] = pwa_estimator( fnonlin, xv, uv, C, D, flag )
%UNTITLED5 Summary of this function goes here
%   flag = 'l' for lipshitz estimator or 'p' for projection

%% selecting estimator
if flag == 'l'
    estimF = @pwa_lipschitz;
else
    estimF = @pwa_projection;
end

%% performing pwa estimation for each state
nx = length(xv);
nu = length(uv);
for i = 1:nx
    fTemp = @(x,u) (1:nx==i)*fnonlin(x,u);
    pwaTemp = estimF(fTemp,xv,uv);
    if i == 1
        post_pwa = pwaTemp;
    else
        for j = 1:numel(post_pwa)
            post_pwa{j}.Bu = [post_pwa{j}.Bu;pwaTemp{j}.Bu];
            post_pwa{j}.A = [post_pwa{j}.A;pwaTemp{j}.A];
            post_pwa{j}.Bw = [post_pwa{j}.Bw;pwaTemp{j}.Bw];
            post_pwa{j}.f = [post_pwa{j}.f;pwaTemp{j}.f];
        end
    end
end
for j = 1:numel(post_pwa)
    post_pwa{j}.B = [post_pwa{j}.A(:,nx+1:end) repmat(post_pwa{j}.Bw,1,nx) post_pwa{j}.f]; % n x (m + n +1)
    post_pwa{j}.A = post_pwa{j}.A(:,1:nx);
    post_pwa{j}.C = C;
    post_pwa{j}.D = D;
    post_pwa{j}.repeatedA = 0;
    post_pwa{j}.repeatedB = 0;
end
% post_pwa = pwa_process(post_pwa);
%     function out = nestedT(varargin)
%         out = fnonlin(varargin{:});
%         out = out(i);
%     end
% end

%% Plot if possible
% if (nx+nu==2)
% pwa =reshape(post_pwa,[],1);
% f=fnonlin{1};
%     figure;
%     xmin = xv{1}(1); xmax = xv{1}(end);
%     umin = uv{1}(1); umax = uv{1}(end);
%     [Xmesh,Ymesh]=meshgrid(linspace(xmin(1),xmax(1),20)',linspace(umin(1),umax(1),20)');
%     Z = f(Xmesh,Ymesh);
%     surf(Xmesh,Ymesh,Z);
%
%     figure;
%     hold all
%     [Xmesh,Ymesh]=meshgrid(linspace(xmin(1),xmax(1),50)',linspace(umin(1),umax(1),50)');
%     Z = f(Xmesh,Ymesh);
%     %         surf(Xmesh,Ymesh,Z);
%     for ii=1:(size(post_pwa,1)*size(post_pwa,2))
%         [Xmesh,Ymesh]=meshgrid(linspace(pwa{ii}.lim(1,1),pwa{ii}.lim(2,1),2)',linspace(pwa{ii}.lim(1,2),pwa{ii}.lim(2,2),2)');
%         Z = pwa{ii}.A*Xmesh + pwa{ii}.Bu*Ymesh + pwa{ii}.f;
%         surf(Xmesh,Ymesh,Z);
%         Z = pwa{ii}.A*Xmesh + pwa{ii}.Bu*Ymesh + pwa{ii}.Bw + pwa{ii}.f;
%         surf(Xmesh,Ymesh,Z);
%         Z = pwa{ii}.A*Xmesh + pwa{ii}.Bu*Ymesh - pwa{ii}.Bw + pwa{ii}.f;
%         surf(Xmesh,Ymesh,Z);
%     end
%
% end