function y = hybridBolusFun(t,t0,a,tMax)
    % https://www.desmos.com/calculator/whvydauxkr
    [~,x0] = mink(abs(t - t0),1);
    yP = real(tMax ^ (-a) .* (t - t0) .^ a .* exp(a .* (1 - (t - t0) / tMax)));
    y = [zeros(x0 - 1,1); yP(x0:end)];
end
