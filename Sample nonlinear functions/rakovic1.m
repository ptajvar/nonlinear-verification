function handle = rakovic1()
A{1} = [1 0.2;0 1];
A{2} = [0.5 0.2;0 1];
B = [0;1];
f{1} = [0;0];
f{2} = [0.5;0];

handle = @twoD_lin;

function xn_next = twoD_lin(x,u)
    mod = (x(1)>1)+1;
    xn_next = A{mod}*x + B*u + f{mod};
end

end