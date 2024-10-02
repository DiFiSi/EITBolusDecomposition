function [Yhat,A,X] = YFun(t,p,Y,n,A)
    if isempty(A)
        X = XFun(t,p,n);
        A = AFun(X,Y);
    else
        X = XScaledFun(p,t);
    end
    Yhat = X * A;
end
