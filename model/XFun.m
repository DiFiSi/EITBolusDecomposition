function X = XFun(t,p,n)
    X = zeros(length(t),n);
    for i = 1:n
        idx = (i - 1) * 3;
        X(:,i) = hybridBolusFun(t,p(idx + 1),p(idx + 2),p(idx + 3));
    end
end
