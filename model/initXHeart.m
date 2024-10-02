function [p0, lb, ub] = initXHeart(y,t,TOA,movWin)
    yLow    = smoothdata(y,'gaussian',movWin);
    [yMTT,~,ytSlopeL,ytSlopeR] = fwhm(yLow,t);
    yLenFourth = (ytSlopeR - ytSlopeL) / 4;

    %% RH portion
    % tMax
    [~, rhxMax] = maxk(yLow,1);
    rhMax   = min(y(rhxMax) * [0.20;0.50;1.10],1);
    rhtMax  = [ytSlopeL;ytSlopeL + yLenFourth;ytSlopeR - yLenFourth];
    
    % t0
    rht0    = [0;TOA;ytSlopeL]; 
    rhtMaxCands  = rhtMax - rht0'; 
    rhtMaxCands  = rhtMaxCands(rhtMaxCands > 0);
    rhtMax  = [max(rhtMaxCands);rhtMax(2) - rht0(2);min(rhtMaxCands)];
    
    % alpha
    rha     = [1.00; rhtMax ./ (0.25 * yMTT / 2); 10.00];
    
    %% LH portion
    % tMax
    [~, lhxMax] = maxk(yLow,1);
    lhMax   = min(y(lhxMax) * [0.20;0.50;1.10],1);
    lhtMax  = [ytSlopeL + yLenFourth;ytSlopeR - yLenFourth;ytSlopeR];
    
    % t0
    lht0    = [ytSlopeL;ytSlopeL + yLenFourth;ytSlopeR - yLenFourth]; 
    lhtMaxCands  = lhtMax - lht0'; 
    lhtMaxCands  = lhtMaxCands(lhtMaxCands > 0);
    lhtMax  = [max(lhtMaxCands);lhtMax(2) - lht0(2);min(lhtMaxCands)];
    
    % alpha
    lha     = [1.00; lhtMax ./ (0.25 * yMTT / 2); 10.00];
    
    % Merge p0, lb, and ub
    p0    = [[rhMax(2);   rht0(2);    rha(2);    rhtMax(2)],...
             [lhMax(2);   lht0(2);    lha(2);    lhtMax(2)]];
    lb    = [[rhMax(1);   rht0(1);    rha(1);    rhtMax(3)],...
             [lhMax(1);   lht0(1);    lha(1);    lhtMax(3)]];
    ub    = [[rhMax(3);   rht0(3);    rha(3);    rhtMax(1)],...
             [lhMax(3);   lht0(3);    lha(3);    lhtMax(1)]];
end
