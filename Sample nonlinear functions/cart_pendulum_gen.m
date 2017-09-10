function handle = cart_pendulum_gen(T)
if nargin < 1
    T = 0.1;
end

handle = @cart_pendulum;

function xn_next = cart_pendulum(x,u)
xn_next = [x(1)+T*x(2);x(2)-T*(sin(x(1))+u*cos(x(1)))];
end

end