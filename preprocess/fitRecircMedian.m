function yR = fitRecircMedian(y,xStart,xEnd,disp)
    if ~exist('disp','var')
        disp = 0;
    end
    
    n = length(y);
    x = (1:n)';
    
    % Cubic interpolation between pulse onset and ending to approximate 
    % baseline drift due to recirculation 
    xBefore = unique([1;xStart]);
    xAfter = unique([xEnd;n]);
    xDrift = [xBefore;xAfter];
    yDrift = [repmat(median(y(1:xStart)),length(xBefore),1);repmat(median(y(xEnd:min(n,xEnd+100))),length(xAfter),1)];
    
    yR = interp1(xDrift,yDrift,x,'pchip');
    
    if disp
        figure;
        plot(x,y); hold on;
        scatter(xDrift,yDrift);
        plot(x,yR);
        plot(x,y - yR + min(yR));
    end
end