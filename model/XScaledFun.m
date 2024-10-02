function X = XScaledFun(p,t)
    n = length(p) / 4;
    X = zeros(length(t),n);
    for i = 1:n
        idx = (i - 1) * 4;
        X(:,i) = p(idx + 1) * hybridBolusFun(t,p(idx + 2),p(idx + 3),p(idx + 4));
    end
end
