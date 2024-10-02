function [imgClean,params,paramsNames,mask,giLocs] = applyCorrection(imgFilt, params, paramsNames, mask, fs)
    tSize = size(imgFilt,3);
    t = (0:tSize - 1)' / fs;

    % Spatial smoothing hyperparameters
    gaussWidth = 9;
    gaussSigma = 2;
    gaussFilt = fspecial('gaussian',gaussWidth,gaussSigma);

    % Smoothen volatile parameters
    intIdxs = find(ismember(paramsNames,{'RFNum','B','TOANum','TOENum'}));
    roundIdxs = find(ismember(paramsNames,{'TOANum','TOENum'}));

    nLdmks = 3;
    for i = 1:nLdmks
        idx = intIdxs(i);
        params{idx} = nanconv(params{idx},gaussFilt,'nanout');
        if ismember(idx,roundIdxs)
            params{idx} = round(params{idx});
        end
    end
    
    % Retrieve relevant params for pixel-wise signal correction
    ldmk = cat(3,params{ismember(paramsNames,{'TOANum','TOENum','M','B'})});
    
    % Initialize params
    imgClean = nan(size(imgFilt));
    fail = zeros(size(mask));
    [rows,cols] = findND(mask);
    nIter = length(rows);
    count = 0;
    for i = 1:nIter
        try
            row = rows(i);
            col = cols(i);
            y = squeeze(imgFilt(row,col,:));

            % Detrend whole pixel signal based on pre-pulse segment
            m = ldmk(row,col,3);
            b = ldmk(row,col,4);
            yTrend = m * t + b;
            y = y - yTrend;

            % Remove recirculation using previously found landmarks
            xStart = ldmk(row,col,1);
            xEnd = ldmk(row,col,2);
            yR = fitRecircMedian(y,xStart,xEnd);
            y = y - yR;

            imgClean(row,col,:) = y;

            % Display iteration number for quality control
            count = count + 1;
            disp("Voxel " + num2str(count) + " / " + num2str(nIter));
        catch
            count = count + 1;
            fail(row,col) = 1;
            disp("FAILED Voxel " + num2str(count) + " / " + num2str(nIter));
            continue
        end
    end

    % Bolus localization
    B2P = params{ismember(paramsNames,'B2P')};
    locsMask = (B2P >= 0.20 * max(B2P(:),[],'omitnan'));
    TOANum = params{ismember(paramsNames,"TOANum")};
    TOANum(~locsMask) = NaN;
    TOENum = params{ismember(paramsNames,"TOENum")};
    TOENum(~locsMask) = NaN;
    
    giLocs = zeros(2,1);
    giLocs(1) = ceil(quantile(TOANum(:),0.05));
    giLocs(2) = floor(quantile(TOENum(:),0.95));

    % Update mask
    mask = mask & ~fail;
    nLdmks = length(paramsNames);
    for i = 1:nLdmks
        params{i}(~mask) = NaN;
    end
end
