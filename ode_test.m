function [y] = test()

x0 = [1 3 5
    -1 -3 -5];

function [yf] = grad(t, y)
    size(y)
    yf = -y;
end

[t, y] = ode45(@grad, [0 10], x0);

plot(t, y(:, 1), t, y(:, 2))

end