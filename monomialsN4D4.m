function b = monomialsN4D4(x)
    r = size(x, 1);
    b = zeros(r, s(4, 4));
    for i = 1:r
        b(i, :) = [1, x(i,1), x(i,2), x(i,3), x(i,4), x(i,1)^2, x(i,1)*x(i,2), x(i,2)^2, x(i,1)*x(i,3), x(i,2)*x(i,3), x(i,3)^2, x(i,1)*x(i,4), x(i,2)*x(i,4), x(i,3)*x(i,4), x(i,4)^2, x(i,1)^3, x(i,1)^2*x(i,2), x(i,1)*x(i,2)^2, x(i,2)^3, x(i,1)^2*x(i,3), x(i,1)*x(i,2)*x(i,3), x(i,2)^2*x(i,3), x(i,1)*x(i,3)^2, x(i,2)*x(i,3)^2, x(i,3)^3, x(i,1)^2*x(i,4), x(i,1)*x(i,2)*x(i,4), x(i,2)^2*x(i,4), x(i,1)*x(i,3)*x(i,4), x(i,2)*x(i,3)*x(i,4), x(i,3)^2*x(i,4), x(i,1)*x(i,4)^2, x(i,2)*x(i,4)^2, x(i,3)*x(i,4)^2, x(i,4)^3, x(i,1)^4, x(i,1)^3*x(i,2), x(i,1)^2*x(i,2)^2, x(i,1)*x(i,2)^3, x(i,2)^4, x(i,1)^3*x(i,3), x(i,1)^2*x(i,2)*x(i,3), x(i,1)*x(i,2)^2*x(i,3), x(i,2)^3*x(i,3), x(i,1)^2*x(i,3)^2, x(i,1)*x(i,2)*x(i,3)^2, x(i,2)^2*x(i,3)^2, x(i,1)*x(i,3)^3, x(i,2)*x(i,3)^3, x(i,3)^4, x(i,1)^3*x(i,4), x(i,1)^2*x(i,2)*x(i,4), x(i,1)*x(i,2)^2*x(i,4), x(i,2)^3*x(i,4), x(i,1)^2*x(i,3)*x(i,4), x(i,1)*x(i,2)*x(i,3)*x(i,4), x(i,2)^2*x(i,3)*x(i,4), x(i,1)*x(i,3)^2*x(i,4), x(i,2)*x(i,3)^2*x(i,4), x(i,3)^3*x(i,4), x(i,1)^2*x(i,4)^2, x(i,1)*x(i,2)*x(i,4)^2, x(i,2)^2*x(i,4)^2, x(i,1)*x(i,3)*x(i,4)^2, x(i,2)*x(i,3)*x(i,4)^2, x(i,3)^2*x(i,4)^2, x(i,1)*x(i,4)^3, x(i,2)*x(i,4)^3, x(i,3)*x(i,4)^3, x(i,4)^4];
    end
end