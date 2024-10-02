function [Aout,p,err] = fitA(Y,X,p0,lb,ub)
    % Decompose lung compartment into two compartments: right and left lung
    % assuming it is well approximated by a swarm of small bivariate
    % Gaussians. After fitting said Gaussians, they are sorted into the 
    % right and left lung compartments and added up.

    % Cost function
    AFun = @normmodel;

    % Set Gaussian swarm initial guesses and bounds
    p0 = p0(1:end-2,2:3);
    p0(4:5,:) = repmat([0.05,0.05],2,1);

    lb = lb(1:end-2,2:3);
    lb(3,:) = [0,0.5];
    lb(4:5,:) = repmat([0,0],2,1);
    lb(6,:) = [0,0];

    ub = ub(1:end-2,2:3);
    ub(3,:) = [0.5,1];
    ub(4:5,:) = repmat([0.1,0.1],2,1);
    ub(6,:) = [0,0];

    order = 5;
    p0 = repmat(p0,1,order);
    lb = repmat(lb,1,order);
    ub = repmat(ub,1,order);
    
    % Cost function
    n = size(p0,2);
    cost = @(p) sum(AFun(X,p,n),2) - Y;

    % Optimize cost function w/ lsqnonlin
    opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective','StepTolerance',10e-3,'Display','iter');
    [p, err, ~] = lsqnonlin(cost,p0(:),lb(:),ub(:),opts);
    
    % Collect results
    Aout = AFun(X,p,n);
end
