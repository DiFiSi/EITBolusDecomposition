function [xStart, xEnd, xMax] = findBolus(imgFilt)
    % Pre-trim general bolus region after finding start and end
    ySum = squeeze(sum(imgFilt,[1,2],'omitnan'));
    tSize = length(ySum);
    
    % Get maximum of global signal
    y = ySum;
    [~,xMax] = findpeaks(y,"NPeaks",1,"SortStr","descend","MinPeakProminence",range(y) * 0.2);

    % Run Heron formula on global signal to estimate pulse onset and ending
    xStart = findStartEndBolusLims(y, [1,xMax]);
    xEnd = findStartEndBolusLims(y, [xMax,tSize]);
end