close all
f = delay_logistics_gen;
uv = cell(0,0);
xv = cell(1,2);
xv{1} = [1.1:.5:2.6];
xv{2} = [1.1:0.5:2.6];

xv2 = cell(1,2);
xv2{1} = [1.1:.3:2.6];
xv2{2} = [1.1:0.3:2.6];

pwa = pwa_estimator( f, xv, uv, 0, 0, 'l' );
pwa2 = pwa_estimator( f, xv2, uv, 0, 0, 'l' );

nt = 100;
x = zeros(2,nt);
x(:,1) = [1.2;1.5];
err = zeros(1,nt);
err2 = zeros(1,nt);
xPWA = x;
xPWA2 = x;
for i = 2:nt
    x(:,i) = f(x(:,i-1));
    [xPWA(:,i),E] = eval_pwa(pwa,x(:,i-1));
    [xPWA2(:,i),E2] = eval_pwa(pwa2,x(:,i-1));
    err(i) = max(abs(E(1,:)));
    err2(i) = max(abs(E2(1,:)));
end
subplot(1,2,1)
plot(x(1,:),'LineWidth',2);
hold on
plot(xPWA(1,:)+err,'LineWidth',2);
plot(xPWA(1,:)-err,'LineWidth',2);
title('(a)')
xlabel('t');
ylabel('x');
for i = 1:length(xv{2})
    plot([0 nt],[xv{2}(i) xv{2}(i)],'k')    
end

subplot(1,2,2)
plot(x(1,:),'LineWidth',2);
hold on
plot(xPWA2(1,:)+err2,'LineWidth',2);
plot(xPWA2(1,:)-err2,'LineWidth',2);
title('(b)')

legend('Delay Logistics Function', 'PWA approximation upper bound', 'PWA approximation lower bound');


for i = 1:length(xv2{2})
    plot([0 nt],[xv2{2}(i) xv2{2}(i)],'k')    
end

xlabel('t');
ylabel('x');


set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 16 10]); %x_width=10cm y_width=15cm

print('delayLogistics.png','-dpng','-r350');
