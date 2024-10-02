%% Set directories
addpath("utils\");
addpath("preprocess\");
addpath("model\");
addpath("other\");

%% Load data
load("data.mat");

img     = data.img;
t       = data.t;
fs      = data.fs;
dLp     = data.dLp;
ctROIs  = data.ctROIs;

%% Pre-processing
% Low-pass filter to isolate bolus pulse
imgFilt = mirroredFilt(img,dLp);

% Bolus localization, landmark identification (onset, ending, maxima), and
% preliminary parameter estimation
[giStart, giEnd, giMax] = findBolus(imgFilt);
[params, paramsNames, mask] = landmarkBolus(imgFilt, giMax, [giStart;giEnd], ctROIs > 0, fs);

% Correcting the signal using identified landamrks (e.g., removal of 
% baseline drift due to recirculation)
[imgCorr,params,paramsNames,mask,giLocs] = applyCorrection(imgFilt, params, paramsNames, mask, fs);

% Get mask indexes for later
[xSize,ySize,~] = size(imgFilt);
[maskRows,maskCols] = find(mask);
maskIdxs = sub2ind([xSize,ySize],maskRows,maskCols);

% Trim entire corrected image between onset and ending of global signal
imgTrim = imgCorr(:,:,giLocs(1):min(giLocs(2) + fs * 3,size(imgCorr,3)));
tTrim = (0:size(imgTrim,3) - 1)' / fs;

%% Plotting results of pre-processing
yGlobal = squeeze(sum(img,[1,2],'omitnan'));
yGlobalFilt = squeeze(sum(imgFilt,[1,2],'omitnan'));
yGlobalCorr = squeeze(sum(imgCorr,[1,2],'omitnan'));

% Processed curve
figure; hold on;
sgtitle("Preprocessing Overview");

plot(t,yGlobal,'LineWidth',1,'Color',[0.75,0.75,0.75]);
plot(t,yGlobalFilt,'LineWidth',1,'Color',[0.5,0.5,0.5]);
plot(t,yGlobalCorr,'LineWidth',2,'Color',[0,0,0]);
plot(t(giLocs(1)),yGlobalCorr(giLocs(1)),...
     'Marker','>','MarkerSize',10,'LineWidth',2,'Color','k');
plot(t(giLocs(2)),yGlobalCorr(giLocs(2)),...
     'Marker','<','MarkerSize',10,'LineWidth',2,'Color','k');
plot(t(giMax),yGlobalCorr(giMax),...
     'Marker','o','MarkerSize',10,'LineWidth',2,'Color','k');
xlabel("Time [s]");
ylabel("Ampltiude [a.u.]");
legend(["Raw Global","Filt. Global","Corrected Global","Onset","End","Peak"]);

% Preliminary parameters
figure;
sgtitle("Prelinimary Metrics");

subplot(3,2,1);
imagesc(params{ismember(paramsNames,'TOANum')} ./ fs);
axis off;
colorbar;
title("Time of Arrival [s]");

subplot(3,2,2);
imagesc(params{ismember(paramsNames,'TOENum')} ./ fs);
axis off;
colorbar;
title("Time of Ending [s]");

subplot(3,2,3);
imagesc(params{ismember(paramsNames,'M')});
axis off;
colorbar;
title("Apnea Drift Slope [a.u. / s]");

subplot(3,2,4);
imagesc(params{ismember(paramsNames,'B')});
axis off;
colorbar;
title("Apnea Drift Intercept [a.u.]");

subplot(3,2,5);
imagesc(params{ismember(paramsNames,'B2P')});
axis off;
colorbar;
title("Baseline-to-peak [a.u.]");

subplot(3,2,6);
imagesc(params{ismember(paramsNames,'RFNum')});
axis off;
colorbar;
title("Recirculation Amplitude [a.u.]");

%% Global NMF 
% General decomposition formulation Y = XA
% Y [time x #pixels]            - data
% X [time x #compartments]      - temporal dynamics
% A [#compartments x #pixels]   - spatial dynamics

% Normalize signal for better optimization computation
imgTrimNorm = imgTrim / max(imgTrim(:),[],'omitnan');
imgRef = max(imgTrim,[],3);
imgRefNorm = imgRef / max(imgTrim(:),[],'omitnan');

% Calculcate initial guesses and bounds for NMF parameters (usually not 
% game-changing, but improves speed at least)
[p0X,lbX,ubX,p0A,lbA,ubA] = initXA(imgTrimNorm,params,paramsNames,ctROIs,giLocs,fs);

% Set data
Y = reshape(imgTrimNorm,xSize * ySize,[])';
Y = Y(:,maskIdxs);
Y(isnan(Y)) = 0;

% Run optimization for 3 compartments (lungs joined into one compartment)
p0 = [p0X(:,1), mean(p0X(:,2:3),2), p0X(:,4)];
lb = [lbX(:,1), mean(lbX(:,2:3),2), lbX(:,4)];
ub = [ubX(:,1), mean(ubX(:,2:3),2), ubX(:,4)];
[Ytmp,Xtmp,Atmp,pTmp] = fitNMF(Y,tTrim,p0(:),lb(:),ub(:));

imgAtmp = imgFy(Atmp', maskIdxs, [xSize,ySize,3]);

%% Plot results of 3-compartment decomposition
% Spatial and temporal dynamics
figure;
sgtitle("Spatial and Temporal Dynamics");

subplot(3,3,1);
imagesc(imgAtmp(:,:,1));
axis off;
colorbar;
title("Spatial Map Pre-Lung");

subplot(3,3,4);
imagesc(imgAtmp(:,:,2));
axis off;
colorbar;
title("Spatial Map Lungs");

subplot(3,3,7);
imagesc(imgAtmp(:,:,3));
axis off;
colorbar;
title("Spatial Map Post-Lung");

subplot(3,3,[2,3,5,6,8,9]);
plot(tTrim,Xtmp,'LineWidth',2);
xlabel("Time [s]");
ylabel("Ampltiude [a.u.]");
legend(["Pre-Lung","Lungs","Post-Lung"]);
title("Temporal Dynamics");
grid on;

%% Spatial decomposition of lung spatial templates A0
% Set normalized coordinates and lung compartment to decompose
normCoords = [maskRows(:) maskCols(:)] ./ max(xSize,ySize);
A = Atmp(2,:)';

% Run decomposition of lung compartment
[Alung,pLung] = fitA(A,normCoords,p0A,lbA,ubA);
pLungMat = reshape(pLung,6,[]);

% Assign fitted Gaussians to new right and left lung compartments
gaussIdx = pLungMat(3,:) < 0.5;
Alung = [sum(Alung(:,gaussIdx),2),sum(Alung(:,~gaussIdx),2)];

% Normalize spatial templates for right heart, right lung, left lung, and
% left heart
A0 = [Atmp(1,:);Alung';Atmp(3,:)];
imgA0 = imgFy(A0', maskIdxs, [xSize,ySize,4]);
A0 = A0 ./ max(A0,[],2,'omitnan');

%% Plot results of lung compartment decomposition
% Separate lung templates
figure;
sgtitle("Lung Compartment Decomposition");

subplot(2,3,[1,2,4,5]);
imagesc(imgAtmp(:,:,2));
axis off;
colorbar;
title("Spatial Map Lungs");

subplot(2,3,3);
imagesc(imgA0(:,:,2));
axis off;
colorbar;
title("Spatial Map Right Lung");

subplot(2,3,6);
imagesc(imgA0(:,:,3));
axis off;
colorbar;
title("Spatial Map Left Lung");

%% Run final NMF optimization for 4 compartments
p0 = [p0A(1,:);p0X];
lb = [lbA(1,:);lbX];
ub = [ubA(1,:);ubX];
[Y,X,A0,p] = fitNMF(Y,tTrim,p0(:),lb(:),ub(:),A0);

% Rescale A0 to estimated compartment amplitudes
pMat = reshape(p,4,[]);
A = pMat(1,:) .* A0';
X = X ./ pMat(1,:);
imgA = imgFy(A, maskIdxs, [xSize,ySize,4]);

% Obtain multi-compartment ROI
[params,paramsNames] = decidePeaksNMF(imgA,params,paramsNames);
nPeaksMap = sum(params{ismember(paramsNames,'NPEAKS')},3);

%% Plot results of 4-compartment decomposition
% Spatial and temporal dynamics
figure;
sgtitle("Spatial and Temporal Dynamics");

subplot(2,4,1);
imagesc(imgA(:,:,1));
axis off;
colorbar;
title("Spatial Map RH");

subplot(2,4,2);
imagesc(imgA(:,:,4));
axis off;
colorbar;
title("Spatial Map LH");

subplot(2,4,5);
imagesc(imgA(:,:,2));
axis off;
colorbar;
title("Spatial Map RL");

subplot(2,4,6);
imagesc(imgA(:,:,3));
axis off;
colorbar;
title("Spatial Map LL");

subplot(2,4,[3,4,7,8]);
plot(tTrim,X,'LineWidth',2);
xlabel("Time [s]");
ylabel("Ampltiude [a.u.]");
legend(["RH","RL","LL","LH"]);
title("Temporal Dynamics");
grid on;

% Multi-compartment ROI
figure;
sgtitle("Number of Compartments");

imagesc(nPeaksMap);
axis off;
colormap(lines(5));
cb = colorbar;
cb.Ticks = [0,1,2,3,4];
cb.TickLabels = ["0","1","2","3","4"];
