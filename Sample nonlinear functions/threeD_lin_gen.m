function handle = threeD_lin_gen(A,B)
if nargin<2
    B = [0;0;1];
    if nargin<1
        A = [1 1 0;0 1 1;0 0 1];
    end
end

handle = @threeD_lin;

function xn_next = threeD_lin(x,u)
    xn_next = A*x + B*u;
end

end