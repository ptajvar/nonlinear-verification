function handle = rakovic2()
A{1} = [1 1;0 1];
A{2} = [0.5 0.2;0 1];
A{3} = [1 -1;0 1];
A{4} = [1 -1;0 1];

B{1} = [1;0.5];
B{2} = [-1;-0.5];
B{3} = [-1;0.5];
B{4} = [1;-0.5];

MODE = [2 3;4 1];

handle = @func;

function xn_next = func(x,u)
    mod = MODE((x(1)>0)+1,(x(2)>0)+1);
    xn_next = A{mod}*x + B{mod}*u;
end

end