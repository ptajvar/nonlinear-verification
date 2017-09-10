function handle = twoD_rotate(A,B)
if nargin<2
    B = [-1;1];
    if nargin<1
        A = [0.8 .1;-.1 .8];
    end
end

handle = @twoD_lin;

function xn_next = twoD_lin(x,u)
    xn_next = A*x + B*u +0.1*rand(2,1);
end

end