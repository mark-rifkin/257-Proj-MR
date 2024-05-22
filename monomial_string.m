function out = monomial_string(n, d)
    out = eval_monomials(sym('x', [1 n]), d);
    out = append("[", strjoin(string(out), ","), "]");
    for j = 1:n
        out = strrep(out, append("x", num2str(j)), append("x(i,",  num2str(j), ")"));
    end
end