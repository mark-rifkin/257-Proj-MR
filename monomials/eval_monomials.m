% Evaluate standard monomials of vector x up to degree d
% Recursion in monomial_powers is very slow. To generate hard coded
% monomials: use monomial_string and write into method like monomialsN4D4
function b = eval_monomials(x, d)
   r = size(x, 1);
   n = size(x, 2);

   % hardcoded cases
   if d == 1
       b = [ones(r, 1), x]; return;
   elseif n == 4
       if d == 4
           b = monomialsN4D4(x); return;
       elseif d == 2
           b = monomialsN4D2(x); return;
       end
   elseif n == 3
       if d == 2
           b = monomialsN3D2(x); return
       end
   end

   

   if isnumeric(x)
      b = ones(r, s(n, d));
   else
      b = sym(ones(r, s(n, d)));
   end

   for row = 1:r
       col = 1;
       for degree = 1:d % iterate over each total degree
           powers = monomial_powers(n, degree);

           for i = 1:size(powers, 1)
               for j = 1:n
                    b(row, col + i) = ...
                        b(row, col+i) * x(row, j)^powers(i, j);
               end
           end

           col = col + size(powers, 1);  
       end
   end
end