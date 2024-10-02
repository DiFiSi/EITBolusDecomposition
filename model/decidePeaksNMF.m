function [params,paramsNames] = decidePeaksNMF(A0,params,paramsNames)
    threshSum = 0.05;
    threshAmp = 0.05;  
    
    % Calculate relative peak importance maps
    maps = abs(A0);
    mapsNorm = maps ./ repmat(sum(maps,3,'omitnan'),1,1,size(maps,3));
    
    % Get NPEAKS map
    NPEAKS = mapsNorm > threshSum & maps > threshAmp;
    params{end+1,1} = NPEAKS;
    paramsNames{end+1,1} = 'NPEAKS';
end