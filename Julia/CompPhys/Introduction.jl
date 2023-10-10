a = [1, 0.1];
b = [-3, -3];
c = [1, 0.1];

function QuadraticEquation(a, b, c)
q = -0.5 .* (b[:] + sign.(b[:]) .* sqrt.(b.^2 - 4 .* a .* c));
x_1 = q./a;
x_2 = c./q;
end

QuadraticEquation(a, b, c)