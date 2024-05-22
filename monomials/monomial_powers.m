% Recursively generate all possible nonnegative integer vectors
% of length n with fixed sum (monomial degree)
% based on https://homepages.laas.fr/henrion/software/momgraph/powvec.m
% but rewritten to use standard monomial ordering
function powers = monomial_powers(n, sum)
    if n > 1
        powers = zeros(1, n);
        if sum > 0
            row = 1;
            for end_power = 0:sum
                r_powers = monomial_powers(n-1, sum-end_power);
                r_l = size(r_powers, 1);
                powers(row:row+r_l-1, end) = repmat(end_power,r_l, 1);
                powers(row:row+r_l-1, 1:end-1) = r_powers;
                
                row = row+r_l; 
            end
        end
    else
        powers = sum;
    end
end
  