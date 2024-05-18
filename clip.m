function x = clip(x, lb, ub)
    x = min(max(x, lb), ub);
end
