function b_eval = monomialsN4D2(x) % input x is a r by 4 matrix of points
    r = size(x, 1);
    b_eval = zeros(r, s(4, 2));
    for i = 1:r
        b_eval(i, :) = [1, x(i, :), ...
        x(i, 1)^2, x(i, 1)*x(i, 2), x(i, 2)^2,  ...
        x(i, 1)*x(i, 3), x(i, 2)*x(i, 3), x(i, 3)^2, ...
        x(i, 1)*x(i, 4), x(i, 2)*x(i, 4), x(i, 3)*x(i, 4), x(i, 4)^2];
    end
end