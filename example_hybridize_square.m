f = @(x,u) [(x(1,:)+x(2,:)).^2;x(1,:)];
nx = 2;
xv = {[-1 1],[-1 1]};
uv = {};

pwa = pwa_estimator( f, xv, uv, 0, 0, 'l' );

x1 = -1:.1:1;
x2 = -1:.05:1;
[X1,X2] = meshgrid(x1,x2);
res = f([X1(:)';X2(:)']);
res = res(1,:);
res = reshape(res,size(X1));
mesh(X1,X2,res);
hold on
res = eval_pwa(pwa,[X1(:)';X2(:)']);
res = res(1,:);
res = reshape(res,size(X1));
mesh(X1,X2,res);
mesh(X1,X2,res+pwa{1}.Bw(1,1));
mesh(X1,X2,res+pwa{1}.Bw(1,2));