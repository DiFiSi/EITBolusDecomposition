function Aout = normmodel(X,p,n)
    Aout = zeros(size(X,1),n);
    for i = 1:n
        idx = (i - 1) * 6;

        A = p(idx + 1);
        muX = p(idx + 2);
        muY = p(idx + 3);
        scaleX = p(idx + 4);
        scaleY = p(idx + 5);
        rho = p(idx + 6);
        
        Xtrans = shiftScaleX(X,[muX,muY],[scaleX,scaleY]);
        Sigma = [1, rho; rho, 1];
        Aout(:,i) =  A * nativenorm(Xtrans, Sigma);
    end
end
