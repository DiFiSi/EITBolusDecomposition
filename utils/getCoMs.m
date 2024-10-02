function CoMs = getCoMs(ROIs,img)
    % This function provides the centers of mass for all non-zero labelled
    % regions of interest in the ROIs image. ROIs should be a 2D matrix
    % with integer elements as pixel-wise labels

    if exist('img','var')
        doWeighted = 1;
    else
        doWeighted = 0;
    end
    
    labels = unique(ROIs(ROIs > 0));
    nLabels = numel(labels);
    
    CoMs = zeros(nLabels,2);
    
    % Center of mass weighted with signal's standard deviation
    if doWeighted
        for i = 1:nLabels
            ROI = bwareafilt(ROIs == labels(i),1,'largest');
            
            prop = regionprops(ROI,img,'WeightedCentroid');
            CoMs(i,:) = [prop.WeightedCentroid(2),prop.WeightedCentroid(1)];
        end
    % Center of mass without weighting
    else
        for i = 1:nLabels
            ROI = bwareafilt(ROIs == labels(i),1,'largest');
            
            prop = regionprops(ROI,'Centroid');
            CoMs(i,:) = [prop.Centroid(2),prop.Centroid(1)];
        end
    end
    
    CoMs = round(CoMs);
end