function [params, paramsNames, mask] = landmarkBolus(imgFilt, giMax, giLocs, mask, fs)
    % Get data dimensions
    xSize = size(imgFilt,1);
    ySize = size(imgFilt,2);
    tSize = size(imgFilt,3);
    t = (0:tSize - 1)' / fs;
    
    % Processing hyperparameters
    movWin = 5 * fs;
    sigExtend = 5 * fs;

    % Get coordinates of valid pixels
    [rows,cols] = findND(mask);
    nIter = length(rows);
    
    % Initialize output params
    TOANum = nan(xSize,ySize);
    TOENum = nan(xSize,ySize);
    M = nan(xSize,ySize);
    B = nan(xSize,ySize);
    B2P = nan(xSize,ySize);
    RFNum = nan(xSize,ySize);
    
    count = 0;
    fail = zeros(size(mask));
    for i = 1:nIter
        try
            row = rows(i);
            col = cols(i);
            y = squeeze(imgFilt(row,col,:));

            % Detrend whole pixel signal based on pre-pulse segment
            [yDetrend,m,b] = detrendBolusRob(y,t,giLocs(1));

            % Smoothen pixel signal
            yLow = smoothdata(yDetrend,'gaussian',movWin);
            
            % Findpeaks in bounded signal
            [pk,~] = maxk(yLow(giLocs(1):giLocs(2)),1);

            % Find xStart and xEnd with Heron triangle approach
            yStart = yLow(giLocs(1):-1:1);
            baseStart = mean(yStart);
            xStart = findStartEndBolusLims([yLow(1) * ones(sigExtend,1);yLow],...
                     [1,giMax + sigExtend],sigExtend + giLocs(1) + [-2,2] * fs);
            xStart = xStart - sigExtend;
            xEnd = findStartEndBolusLims([yLow;yLow(end) * ones(sigExtend,1)],...
                   [giMax,tSize + sigExtend],giLocs(2) + [-2,2] * fs);
            
            % Get B2P (baseline-to-peak) & maximum recirculation value
            b2p = pk - baseStart;
            RF = max(0,yDetrend(xEnd) - baseStart);

            % Assign vales to params vectors
            TOANum(row,col) = xStart;
            TOENum(row,col) = xEnd;
            M(row,col) = m;
            B(row,col) = b;
            B2P(row,col) = b2p;
            RFNum(row,col) = RF;
            
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
    
    % Organize params
    params{1,1} = TOANum;
    paramsNames{1,1} = 'TOANum';
    params{2,1} = TOENum;
    paramsNames{2,1} = 'TOENum';
    params{3,1} = M;
    paramsNames{3,1} = 'M';
    params{4,1} = B;
    paramsNames{4,1} = 'B';
    params{5,1} = B2P;
    paramsNames{5,1} = 'B2P';
    params{6,1} = RFNum;
    paramsNames{6,1} = 'RFNum';

    mask = mask .* ~fail;
end
