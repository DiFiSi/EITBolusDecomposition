function p = nativenorm(X, Sigma) 
    p = mvnpdf(X, [0, 0], Sigma);
    p = p ./ max(p);
end
