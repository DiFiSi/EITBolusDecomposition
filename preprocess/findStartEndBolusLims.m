function xBMax = findStartEndBolusLims(y, verts, lims)
    % Function to compute sides of triangle
    computeTriangleDistances = @(ax, ay, bx, by, cx, cy) [sqrt((bx-cx).^2+(by-cy).^2),...
                                                          sqrt((ax-cx).^2+(ay-cy).^2),...
                                                          sqrt((ax-bx).^2+(ay-by).^2)];
                    
    % Function to compute triangle area with Heron's formula                                                  
    heronFormula = @(a,b,c) sqrt((a+b+c)./2 .*(((a+b+c)./2)-a).*(((a+b+c)./2)-b).*(((a+b+c)./2)-c));
    
    % Set the leftmost (A) and rightmost (C) vertices of triangle, between
    % which the point B will move
    xA = verts(1); % x-coordinate of vertex A
    xC = verts(2); % x-coordinate of vertex C

    % Number of positions of point B
    nPoints = xC - xA + 1;
    
    % Coordinates of all possible triangle vertices
    xB = (verts(1):verts(2))';
    yB = y(xB);
    xA = repmat(xA,nPoints,1);
    yA = y(xA);
    xC = repmat(xC,nPoints,1);
    yC = y(xC);
    
    % Areas of all possible triangle configurations
    dists = computeTriangleDistances(xA, yA, xB, yB, xC, yC);
    areas = heronFormula(dists(:,1), dists(:,2), dists(:,3)) ./ (xC - xA);
    
    % Choose onset as the x-coordinate of point B which maximizes the
    % triangle area
    if exist('lims','var')
        remIdxs         = xB < lims(1) | xB > lims(2);
        xB(remIdxs)     = [];
        areas(remIdxs)  = [];
    end
    [~, xBMaxI] = max(areas);
    xBMax = xB(xBMaxI);
end

