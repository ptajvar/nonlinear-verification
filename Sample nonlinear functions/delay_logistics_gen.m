function handle = delay_logistics_gen(A,B)
if nargin == 0
    A = 10;
    B = 5;
end

handle = @delay_logistics;

function xn_next = delay_logistics(x,u)
xn_next = [A.*x(1,:)./(1+B.*x(2,:));x(1,:)];
end

end