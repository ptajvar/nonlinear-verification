function [ controlled_sys, controller ] = Add_Controller( pwa_in, control )
%ADD_CONTROLLER inputs a PWA system and a control law in form of
%optimization output and outputs an autonomous PWA system that is the
%result of application of the control law to the main system
%   Detailed explanation goes here

controlled_sys = pwa_in; % autonomous system after applying the controller
S = size(pwa_in);
nx = pwa_in{1}.nx;
nu = pwa_in{1}.nu;
Cgrid = zeros(1,nx+nu);
for i = 1:nu
    Cgrid(i) = length(pwa_in{1}.grid{nx+i});
end
for i = 1:nx
    Cgrid(nu+i) = length(pwa_in{1}.grid{i});
end
controller = cell(Cgrid-1); % the controller pwa system, maps a state to the control action
controller{1}.nu = nx;
controller{1}.nx = nu;
controller{1}.map = [pwa_in{1}.map(nx+1:end) pwa_in{1}.map(1:nx)];
controller{1}.grid = [pwa_in{1}.grid(nx+1:end) pwa_in{1}.grid(1:nx)];
for i = 1:numel(pwa_in)
    controller{i}.f = zeros(nu,1);
    controller{i}.Bw = zeros(nu,2);
    controller{i}.C = zeros(0,nu);
    controller{i}.D = zeros(0,0);    
    [sub{1:nx}] = ind2sub(S,i);
    low = zeros(1,nx);
    high = zeros(1,nx);
    for j = 1:nx
        low(j) = findValue(['LOWER ' num2str(j) ' ' num2str(sub{j})]);
        high(j) = findValue(['HIGHER ' num2str(j) ' ' num2str(sub{j})]);
        for k = 1:nu
            Ac=findValue(['UG_' num2str(i) '_' num2str(k) ' ' num2str(j)])...
                *2/(high(j)-low(j));
            controlled_sys{i}.A(:,j) = controlled_sys{i}.A(:,j)+controlled_sys{i}.Bu(:,k)*Ac;
            controller{i}.Bu(k,j) = Ac;
            controlled_sys{i}.f = controlled_sys{i}.f-controlled_sys{i}.Bu(:,k)*Ac*((high(j)+low(j)))/2;
            controller{i}.f(k) = controller{i}.f(k)-Ac*((high(j)+low(j)))/2;
        end
    end
    for k = 1:nu
        controlled_sys{i}.f = controlled_sys{i}.f +controlled_sys{i}.Bu(:,k)*...
            findValue(['fc_' num2str(i) ' ' num2str(k)]);
        controller{i}.f(k) = controller{i}.f(k) + findValue(['fc_' num2str(i) ' ' num2str(k)]);
    end
    controlled_sys{i}.Bu = zeros(nx,0);
    controller{i}.A = zeros(nu);
end
controlled_sys{1}.nu = 0;

    function value = findValue(name)
        SS = size(control);
        for it = 1:SS(1)
            if strcmp(control{it,1},name)
                value = control{it,2};
                return;
            end
        end
        error(['Variable: ' name ' not found.']);
    end
end

