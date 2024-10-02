function [Yout,Xout,Aout,p,err] = fitNMF(Y,t,p0,lb,ub,A)
    if ~exist("A","var")
        A = [];
        n = numel(p0) / 3;
    else
        n = numel(p0) / 4;
    end
    
    % Cost function
    cost = @(p) YFun(t,p,Y,n,A) - Y;
    
    % Optimize cost function w/ lsqnonlin
    opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective','StepTolerance',10e-3,'Display','iter');
    [p, err, ~] = lsqnonlin(cost,p0,lb,ub,opts);
    
    % Collect results
    [Yout,Aout,Xout] = YFun(t,p,Y,n,A);
end
