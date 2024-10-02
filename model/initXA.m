function [p0X,lbX,ubX,p0A,lbA,ubA] = initXA(imgTrim,params,paramsNames,ctROIs,giLocs,fs)
    xSize = size(imgTrim,1);
    ySize = size(imgTrim,2);
    tSize = size(imgTrim,3);
    t = (0:tSize - 1)' / fs;
    movWin = 5 * fs;

    % Split ROIs into subROIs (top / bottom, left / right) from CT
    % segmentation and use centers of mass as reference pixels
    refPixs = zeros(4,2);
    refPixs(1,:) = getCoMs(ctROIs == 1);
    refPixs(2,:) = getCoMs(ctROIs == 5);
    refPixs(3,:) = getCoMs(ctROIs == 6);
    refPixs(4,:) = getCoMs(ctROIs == 2);

    %% Collect center of mass signals
    rhY     = squeeze(imgTrim(refPixs(1,1),refPixs(1,2),:)); % right heart
    lhY     = squeeze(imgTrim(refPixs(4,1),refPixs(4,2),:)); % left heart
    rlY     = squeeze(imgTrim(refPixs(2,1),refPixs(2,2),:)); % right lung
    llY     = squeeze(imgTrim(refPixs(3,1),refPixs(3,2),:)); % left lung 

    %% Define model function
    modelFun = @XScaledFun;

    %% Overall strategy
    % Start by using initXHeart or initiXLung to initialize parameters with
    % values read directly from signal. 
    % Then fit simple gamma-variate to refine initial guesses for bigger 
    % problem. If fails, simply use values read directly from data.
    %
    % Start by analysing right heart values and use them to bound following
    % compartment values under the assumption that the right heart peak is
    % the earliest, fastest and sharpest peak.

    %% RH signal landmarks
    rhYFit = rhY;
    rhTOA = (params{ismember(paramsNames,"TOANum")}(refPixs(1,1),refPixs(1,2)) - giLocs(1));
    [rhp0, rhlb, rhub] = initXHeart(rhYFit,t,rhTOA / fs,movWin);

    try
        % Cost function
        cost = @(p) modelFun(p(:,1),t) + modelFun(p(:,2),t) - rhY;

        % Optimize cost function w/ lsqnonlin
        opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective','StepTolerance',10e-3,'Display','off');
        [rhp, ~, res] = lsqnonlin(cost,rhp0,rhlb,rhub,opts);
        idxRes = rhY >= 0;

        NE = (sum(rhY(idxRes)) - sum(abs(res(idxRes)))) / sum(rhY(idxRes));
        if NE < 0.70
            a + b;
        end

        rh0X    = rhp(:,1);
    catch
        rh0X    = rhp0(:,1);
    end

    rhlbX   = rhlb(:,1);
    rhubX   = rhub(:,1);

    %% LH signal landmarks
    lhYFit = lhY;
    lhTOA = (params{ismember(paramsNames,"TOANum")}(refPixs(2,1),refPixs(2,2)) - giLocs(1));
    [lhp0, lhlb, lhub] = initXHeart(lhYFit,t,lhTOA / fs,movWin);

    % Fit
    try
        % Cost function
        cost = @(p) modelFun(p(:,1),t) + modelFun(p(:,2),t) - lhY;

        % Optimize cost function w/ lsqnonlin
        opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective','StepTolerance',10e-3,'Display','off');
        [lhp, ~, res] = lsqnonlin(cost,lhp0,lhlb,lhub,opts);
        idxRes = lhY >= 0;

        NE = (sum(lhY(idxRes)) - sum(abs(res(idxRes)))) / sum(lhY(idxRes));
        if NE < 0.70
            a + b;
        end

        lh0X    = lhp(:,2);
    catch
        lh0X    = lhp0(:,2);
    end

    lhlbX   = lhlb(:,2);
    lhubX   = lhub(:,2);

    %% T0 lower bound for lungs from right heart
    yLow = smoothdata(rhY,'gaussian',movWin);
    tSlopeRH = t(find(yLow >= 0.10 * max(yLow),1,'first'));

    %% RL signal landmarks
    rlYFit = rlY;
    rlTOA = (params{ismember(paramsNames,"TOANum")}(refPixs(3,1),refPixs(3,2)) - giLocs(1));
    [rlp0, rllb, rlub] = initXLung(rlYFit,t,tSlopeRH,rlTOA / fs,movWin,fs);
    
    try
        % Cost function
        cost = @(p) modelFun(p,t) - rlY;

        % Optimize cost function w/ lsqnonlin
        opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective','StepTolerance',10e-3,'Display','off');
        [rlp, ~, res] = lsqnonlin(cost,rlp0,rllb,rlub,opts);
        idxRes = rlY >= 0;

        NE = (sum(rlY(idxRes)) - sum(abs(res(idxRes)))) / sum(rlY(idxRes));
        if NE < 0.70
            a + b;
        end

        rl0X    = rlp;
    catch
        rl0X    = rlp0;
    end

    rllbX   = rllb;
    rlubX   = rlub;

    %% LL signal landmarks
    llYFit = llY;
    llTOA = (params{ismember(paramsNames,"TOANum")}(refPixs(4,1),refPixs(4,2)) - giLocs(1));
    [llp0, lllb, llub] = initXLung(llYFit,t,tSlopeRH,llTOA / fs,movWin,fs);

    try
        % Cost function
        cost = @(p) modelFun(p,t) - llY;

        % Optimize cost function w/ lsqnonlin
        opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective','StepTolerance',10e-3,'Display','off');
        [llp, ~, res] = lsqnonlin(cost,llp0,lllb,llub,opts);
        idxRes = llY >= 0;

        NE = (sum(llY(idxRes)) - sum(abs(res(idxRes)))) / sum(llY(idxRes));
        if NE < 0.8
            a + b;
        end

        ll0X    = llp;
    catch
        ll0X    = llp0;
    end

    lllbX   = lllb;
    llubX   = llub;
    
    %% Build bounds X
    p0X     = [rh0X(2:end),    rl0X(2:end),   ll0X(2:end),   lh0X(2:end)];
    lbX     = [rhlbX(2:end),   rllbX(2:end),  lllbX(2:end),  lhlbX(2:end)];
    ubX     = [rhubX(2:end),   rlubX(2:end),  llubX(2:end),  lhubX(2:end)];

    %% Build bounds A
    refPixsNorm = refPixs ./ [xSize,ySize];

    % Assign bounds 
    % Right heart
    rh0A    = [rh0X(1);         refPixsNorm(1,1);   refPixsNorm(1,2);   1;  1;  0;      0;      0   ];
    rhlbA   = [rhlbX(1) * 0.10; 0;                  0;                  0;  0;  -1;     -10;    -10 ];
    rhubA   = [1;               1;                  1;                  1;  1;  1;      10;     10  ];
    
    % Left heart
    lh0A    = [lh0X(1);         refPixsNorm(2,1);   refPixsNorm(2,2);   1;  1;  0;      0;      0   ];
    lhlbA   = [lhlbX(1) * 0.10; 0;                  0;                  0;  0;  -1;     -10;    -10 ];
    lhubA   = [1;               1;                  1;                  1;  1;  1;      10;     10  ];
    
    % Right Lung
    rl0A    = [rl0X(1);         refPixsNorm(3,1);   refPixsNorm(3,2);   1;  1;  0;      0;      0   ];
    rllbA   = [rllbX(1) * 0.10; 0;                  0;                  0;  0;  -1;     -10;    -10 ];
    rlubA   = [1;               1;                  1;                  1;  1;  1;      10;     10  ];
    
    % Left Lung
    ll0A    = [ll0X(1);         refPixsNorm(4,1);   refPixsNorm(4,2);   1;  1;  0;      0;      0   ];
    lllbA   = [lllbX(1) * 0.10; 0;                  0;                  0;  0;  -1;     -10;    -10 ];
    llubA   = [1;               1;                  1;                  1;  1;  1;      10;     10  ];

    % Join A
    p0A     = [rh0A,    rl0A,   ll0A,   lh0A    ];
    lbA     = [rhlbA,   rllbA,  lllbA,  lhlbA   ];
    ubA     = [rhubA,   rlubA,  llubA,  lhubA   ];

    %% Sanity check
    tmp = p0X < lbX;
    p0X(tmp) = lbX(tmp); 

    tmp = p0X > ubX;
    p0X(tmp) = ubX(tmp); 

    tmp = p0A < lbA;
    p0A(tmp) = lbA(tmp); 

    tmp = p0A > ubA;
    p0A(tmp) = lbA(tmp); 

end
