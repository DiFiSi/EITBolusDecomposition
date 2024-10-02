function [p0, lb, ub] = initXLung(y,t,tSlopeRH,TOA,movWin,fs)
    yLow    = smoothdata(y,'gaussian',movWin);
    [MTT,~,tSlopeL,~] = fwhm(y,t);

    % t0
    t0    = sort([tSlopeRH;TOA;tSlopeL]);
    
    % tMax
    [~, xMax] = maxk(yLow,1);
    yMax   = min(y(xMax) * [0.80;1.00;1.20],1);
    tMax  = xMax / fs;
    tMax  = tMax - t0;
    
    % alpha
    a = [1.00; tMax ./ (0.25 * MTT); 10.00];
    
    % Initial bounds
    p0    = [yMax(2);   t0(2);    a(2);    tMax(2)];
    lb    = [yMax(1);   t0(1);    a(1);    tMax(3)];
    ub    = [yMax(3);   t0(3);    a(2);    tMax(1)];
end
