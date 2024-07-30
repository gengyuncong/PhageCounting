close all
clear all
clc

%Data loading and processing
folderPath = [''];
fileNames = {'d20230216_ODbasedMOIresponsesurface_run1_24-Feb-2023.mat',...
    'd20230216_ODbasedMOIresponsesurface_run2_24-Feb-2023.mat'};

data = struct();
for i_file = 1:numel(fileNames)
    tmpData = load([folderPath fileNames{i_file}]);
    tmpData = tmpData.dataSaved;
    fields = fieldnames(tmpData);

    multiple = (i_file-1)*6;
    for i_series = 1:numel(tmpData)
        data(i_series+multiple).run = i_file;

        for i_field = 1:numel(fields)
            data(i_series+multiple).(fields{i_field}) = tmpData(i_series).(fields{i_field});
        end
    end
end

%Remove irrelevant fields for now
data = rmfield(data, {'hill_xScale', 'hill_qMax', 'hill_K1', 'hill_H1', 'hill_K2', 'hill_H2',...
    'weitzdecay_xScale', 'weitzdecay_Q1', 'weitzdecay_Q2', 'weitzdecay_kRate', 'weitzdecay_MOIstar'});

%Convert growth rate (base 2, related to doubling time) to exponential growth rate (base e)
for i_exp = 1:numel(data)
    data(i_exp).growthRateInfected = data(i_exp).growthRateInfected * log(2);
end

%Instructions for plotting
nRepeats = numel(fileNames);
nGrowthRateValues = numel(data)/nRepeats;
runColors = [220 60 60; 60 60 220]/255;
Mspace = logspace(-5, 4, 200);
nspace = 0:1:5;
runMarkers = ['o' 'v'];
legendLocation = {'SouthEast' 'SouthWest'}; %for OD and growth rate, respectively
for i_series = 1:numel(data)
    data(i_series).color = runColors(ceil(i_series/nGrowthRateValues), :);
    data(i_series).marker = runMarkers(ceil(i_series/nGrowthRateValues));
end

CFUconversion = 1e8; %CFU/mL at standard OD = 1.

disp('Preparations are complete.');

%% INDIVIDUAL MOI RESPONSE CURVES BEFORE ADSORPTION CORRECTION
clc
close all

%Preparations for the combined dataset
combinedData = struct();
for i_series = 1:nGrowthRateValues
    combinedData(i_series).odInfectedVector = [];
    combinedData(i_series).growthRateInfectedVector = [];
    combinedData(i_series).odInfected = [];
    combinedData(i_series).growthRateInfected = [];
    combinedData(i_series).MOIseries = [];
    combinedData(i_series).fractionLysogeny = [];
end

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_series = 1:numel(data)
    subplot(2, 3, mod(i_series-1, nGrowthRateValues)+1);
    hold on;

    currentRun = data(i_series).run;

    xToPlot = data(i_series).MOIseries;
    yToPlot = data(i_series).fractionLysogeny(1, :);
    yErrorToPlot = data(i_series).fractionLysogeny(2, :);

    %Visualization
    graphData(i_series) = errorbar(xToPlot, yToPlot, yErrorToPlot,...
        'Marker', 'o', 'MarkerFaceColor', runColors(currentRun, :), 'MarkerEdgeColor', 'None', 'MarkerSize', 4,...
        'LineStyle', 'None', 'Color', runColors(currentRun, :), 'LineWidth', 1, 'CapSize', 6,...
        'DisplayName', ['Infection at OD = ' num2str(data(i_series).odInfected, '%.f')]);

    %For later
    combinedData(mod(i_series-1, nGrowthRateValues)+1).odInfectedVector = [combinedData(mod(i_series-1, 6)+1).odInfectedVector data(i_series).odInfected];
    combinedData(mod(i_series-1, nGrowthRateValues)+1).growthRateInfectedVector = [combinedData(mod(i_series-1, 6)+1).growthRateInfectedVector data(i_series).growthRateInfected(1)];
    combinedData(mod(i_series-1, nGrowthRateValues)+1).MOIseries = [combinedData(mod(i_series-1, 6)+1).MOIseries xToPlot];
    combinedData(mod(i_series-1, nGrowthRateValues)+1).fractionLysogeny = [combinedData(mod(i_series-1, 6)+1).fractionLysogeny [yToPlot; yErrorToPlot]];
end

%Meta-information, and adjusting plot properties
for i_series = 1:numel(combinedData)
    combinedData(i_series).odInfected = [mean([combinedData(i_series).odInfectedVector]); std([combinedData(i_series).odInfectedVector])/sqrt(numel([combinedData(i_series).odInfectedVector]))];
    combinedData(i_series).growthRateInfected = [mean([combinedData(i_series).growthRateInfectedVector]); std([combinedData(i_series).growthRateInfectedVector])/sqrt(numel([combinedData(i_series).growthRateInfectedVector]))];

    subplot(2, 3, i_series);
    if i_series == 1
        legend([graphData(1) graphData(7)], 'Run #1' , 'Run #2', 'Location', 'NorthWest', 'FontSize', 10);
    end
    axis([5e-4 5e3 1e-7 5]);
    xlabel('Bulk MOI (Actual)');
    ylabel('Fraction of lysogeny');
    set(gca, 'XScale', 'Log', 'YScale', 'Log');
    title(['OD = ' num2str(combinedData(i_series).odInfected(1), '%.2f') ' \pm ' num2str(combinedData(i_series).odInfected(2), '%.2f'),...
        ', at ' num2str(combinedData(i_series).growthRateInfected(1), '%.2f') ' \pm ' num2str(combinedData(i_series).growthRateInfected(2), '%.2f') ' /hr'],...
        'FontSize', 11);
    grid on; box on;
end

%% INDIVIDUAL MOI RESPONSE CURVES AFTER CORRECTION FOR ADSORPTION EFFICIENCY À LA MOLDOVAN
clc
close all

%Preparations for the combined dataset
for i_series = 1:numel(combinedData)
    combinedData(i_series).moldovanMOIseries = [];
end

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

%Preparations for Moldovan calculations
CFUconversion = 1e8; %CFU/mL at standard OD = 1.

nBVector = logspace(6, 9.5, 100);
temp = 30; %Celsius
time = 30*60; %seconds of adsorption before dilution
k = 1.6e-11; %For 30C
ksingle = 1.9e-3; %For 30C
kdouble = 7.8e-4; %For 30C

for i_series = 1:numel(data)
    subplot(2, 3, mod(i_series-1, 6)+1);
    hold on;

    currentRun = data(i_series).run;

    %Correcting for the adsorption efficiency
    nB = data(i_series).odInfected * CFUconversion;
    tau1 = 1./(1/2 * ((k.*nB + ksingle + kdouble)...
        + sqrt((k.*nB + ksingle + kdouble).^2 - 4*k.*kdouble.*nB)));
    tau2 = 1./(1/2 * ((k.*nB + ksingle + kdouble)...
        - sqrt((k.*nB + ksingle + kdouble).^2 - 4*k.*kdouble.*nB)));
    fractionFree = 1 ./ (tau2 - tau1) .* ((1/kdouble - tau1) .* exp(-time./tau1)...
        - (1/kdouble - tau2) .* exp(-time./tau2));
    fractionAdsorbed = 1 - fractionFree;

    data(i_series).MOImoldovan = data(i_series).MOIseries * fractionAdsorbed;

    xToPlot = data(i_series).MOImoldovan;
    yToPlot = data(i_series).fractionLysogeny(1, :);
    yErrorToPlot = data(i_series).fractionLysogeny(2, :);

    %Visualization
    graphData(i_series) = errorbar(xToPlot, yToPlot, yErrorToPlot,...
        'Marker', 'o', 'MarkerFaceColor', runColors(currentRun, :), 'MarkerEdgeColor', 'None', 'MarkerSize', 4,...
        'LineStyle', 'None', 'Color', runColors(currentRun, :), 'LineWidth', 1, 'CapSize', 6,...
        'DisplayName', ['Infection at OD = ' num2str(data(i_series).odInfected, '%.f')]);

    %For later
    combinedData(mod(i_series-1, 6)+1).moldovanMOIseries = [combinedData(mod(i_series-1, 6)+1).moldovanMOIseries data(i_series).MOImoldovan];
end

%Meta-information, and adjusting plot properties
for i_series = 1:numel(combinedData)
    subplot(2, 3, i_series);
    if i_series == 1
        legend([graphData(1) graphData(7)], 'Run #1' , 'Run #2', 'Location', 'NorthWest', 'FontSize', 10);
    end
    axis([5e-4 5e3 1e-7 5]);
    xlabel('Bulk MOI (Moldovan-corrected)');
    ylabel('Fraction of lysogeny');
    set(gca, 'XScale', 'Log', 'YScale', 'Log');
    title(['OD = ' num2str(combinedData(i_series).odInfected(1), '%.2f') ' \pm ' num2str(combinedData(i_series).odInfected(2), '%.2f'),...
        ', at ' num2str(combinedData(i_series).growthRateInfected(1), '%.2f') ' \pm ' num2str(combinedData(i_series).growthRateInfected(2), '%.2f') ' /hr'],...
        'FontSize', 11);
    grid on; box on;
end

%% FITTING INDIVIDUAL MOI RESPONSE CURVES - Long computational time
clc
close all

%Preparations for the combined dataset
for i_series = 1:numel(combinedData)
    combinedData(i_series).separateAVector = [];
    combinedData(i_series).separateQ2Vector = [];
    combinedData(i_series).separatePhiVector = [];
    combinedData(i_series).separateKVector = [];
end

%Preparations for fitting
mSpace = logspace(-4, 4, 50);
fitting.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [0 0 0 0], 'Upper', [5 5 1 Inf], 'Startpoint', [0.5 0.1 0.01 0]);
fitting.model = fittype('func_ExpDecaylysogenization3_MOIstar2(a, q2, phi, k, Marray)',...
    'Options', fitting.options, 'Independent', 'Marray', 'Coefficient', {'a', 'q2', 'phi', 'k'});

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_series = 1:numel(data)
    subplot(2, 3, mod(i_series-1, nGrowthRateValues)+1);
    hold on;

    %Data specification
    currentRun = data(i_series).run;
    xToPlot = data(i_series).MOIseries;
    yToPlot = data(i_series).fractionLysogeny(1, :);
    yErrorToPlot = data(i_series).fractionLysogeny(2, :);

    %Fitting
    [xToFit, yToFit] = prepareCurveData(xToPlot, log(yToPlot));
    [param, gof] = fit(xToFit, yToFit, fitting.model);
    fSpace = exp(func_ExpDecaylysogenization3_MOIstar2(param.a, param.q2, param.phi, param.k, mSpace));

    %Visualization
    graphFit(i_series) = plot(mSpace, fSpace,...
        'Marker', 'None', 'LineStyle', '-', 'Color', runColors(currentRun, :), 'LineWidth', 1);

    graphData(i_series) = errorbar(xToPlot, yToPlot, yErrorToPlot,...
        'Marker', runMarkers(currentRun), 'MarkerFaceColor', runColors(currentRun, :), 'MarkerEdgeColor', 'None', 'MarkerSize', 4,...
        'LineStyle', 'None', 'Color', runColors(currentRun, :), 'LineWidth', 1, 'CapSize', 6,...
        'DisplayName', ['Infection at OD = ' num2str(data(i_series).odInfected, '%.f')]);

    %For later
    data(i_series).PQRmodel_xScale = param.a;
    data(i_series).PQRmodel_Q2 = param.q2;
    data(i_series).PQRmodel_phi = param.phi;
    data(i_series).PQRmodel_Q1 = param.q2 * param.phi;
    data(i_series).PQRmodel_k = param.k;

    combinedData(mod(i_series-1, nGrowthRateValues)+1).separateAVector = [combinedData(mod(i_series-1, nGrowthRateValues)+1).separateAVector param.a];
    combinedData(mod(i_series-1, nGrowthRateValues)+1).separateQ2Vector = [combinedData(mod(i_series-1, nGrowthRateValues)+1).separateQ2Vector param.q2];
    combinedData(mod(i_series-1, nGrowthRateValues)+1).separatePhiVector = [combinedData(mod(i_series-1, nGrowthRateValues)+1).separatePhiVector param.phi];
    combinedData(mod(i_series-1, nGrowthRateValues)+1).separateKVector = [combinedData(mod(i_series-1, nGrowthRateValues)+1).separateKVector param.k];
end

%Meta-information, and adjusting plot properties
for i_series = 1:numel(combinedData)
    subplot(2, 3, i_series);
    if i_series == 1
        legend([graphData(1) graphData(7)], 'Run #1' , 'Run #2', 'Location', 'NorthWest', 'FontSize', 10);
    end
    axis([5e-4 5e3 1e-7 5]);
    xlabel('Bulk MOI (Actual)');
    ylabel('Fraction of lysogeny');
    set(gca, 'XScale', 'Log', 'YScale', 'Log');
    title(['OD = ' num2str(combinedData(i_series).odInfected(1), '%.2f') ' \pm ' num2str(combinedData(i_series).odInfected(2), '%.2f'),...
        ', at ' num2str(combinedData(i_series).growthRateInfected(1), '%.2f') ' \pm ' num2str(combinedData(i_series).growthRateInfected(2), '%.2f') ' /hr'],...
        'FontSize', 11);
    grid on; box on;
end

sgtitle('Fitting separately to each experimental replicate');

%% FITTING POOLED MOI RESPONSE CURVES AT THE SAME GROWTH RATE - Long computational time
clc
close all

%Preparations for fitting
mSpace = logspace(-4, 4, 50);
fitting.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [0 0 0 0], 'Startpoint', [0.5 0.1 0.01 0]);
fitting.model = fittype('func_ExpDecaylysogenization3_MOIstar2(a, q2, phi, k, Marray)',...
    'Options', fitting.options, 'Independent', 'Marray', 'Coefficient', {'a', 'q2', 'phi', 'k'});

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

mSpace = logspace(-4, 3.5, 50);
for i_series = 1:numel(combinedData)
    subplot(2, 3, i_series);
    hold on;

    %Data specification
    xToPlot = [data(i_series).MOIseries data(i_series+nGrowthRateValues).MOIseries];
    yToPlot = [data(i_series).fractionLysogeny(1, :) data(i_series+nGrowthRateValues).fractionLysogeny(1, :)];
    yErrorToPlot = [data(i_series).fractionLysogeny(2, :) data(i_series+nGrowthRateValues).fractionLysogeny(2, :)];

    %Fitting
    [xToFit, yToFit] = prepareCurveData(xToPlot, log(yToPlot));
    if i_series == 6 %In stationary phase, constrain for no reduction at high MOI
        [param, gof] = fit(xToFit, yToFit, fitting.model, 'Upper', [5 5 1 0]);
    else
        [param, gof] = fit(xToFit, yToFit, fitting.model, 'Upper', [5 5 1 Inf]);
    end
    fSpace = exp(func_ExpDecaylysogenization3_MOIstar2(param.a, param.q2, param.phi, param.k, mSpace));

    %Visualization
    graphFit(i_series) = plot(mSpace.*param.a, fSpace,...
        'Marker', 'None', 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.2);

    graphData1(i_series) = errorbar(xToPlot(1:6).*param.a, yToPlot(1:6), yErrorToPlot(1:6),...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'MarkerSize', 4,...
        'LineStyle', 'None', 'Color', 'k', 'LineWidth', 1, 'CapSize', 6);

    graphData2(i_series) = errorbar(xToPlot(7:12).*param.a, yToPlot(7:12), yErrorToPlot(7:12),...
        'Marker', 'v', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'MarkerSize', 4,...
        'LineStyle', 'None', 'Color', 'k', 'LineWidth', 1, 'CapSize', 6);

    %Adjusting plot properties
    axis square;
    axis([1e-4 3e2 1e-7 5]);
    xlabel('Bulk MOI');
    ylabel('Fraction of lysogeny');
    set(gca, 'XScale', 'Log', 'YScale', 'Log');
    title(['OD = ' num2str(combinedData(i_series).odInfected(1), '%.2f') ' \pm ' num2str(combinedData(i_series).odInfected(2), '%.2f'),...
        ', at ' num2str(combinedData(i_series).growthRateInfected(1), '%.2f') ' \pm ' num2str(combinedData(i_series).growthRateInfected(2), '%.2f') ' /hr'],...
        'FontSize', 11);
    grid on; box on;

    %For later
%     combinedData(i_series).separateAVector = [combinedData(i_series).separateAVector param.a];
%     combinedData(i_series).separateQ2Vector = [combinedData(i_series).separateQ2Vector param.q2];
%     combinedData(i_series).separatePhiVector = [combinedData(i_series).separatePhiVector param.phi];
%     combinedData(i_series).separateKVector = [combinedData(i_series).separateKVector param.k];
end

sgtitle('Fitting to the pooled data from both experimental replicates');

%% SINGLE-CELL RESPONSE CURVE AT EACH GROWTH RATE
clc
close all

%Preparations
nSpace = 0:1:20;
func_singlecellProb_MOIstar2 = @(q2, phi, k, n) 0 .* (n == 0)...
    + (phi * q2) .* (n > 0 & n < 2)...
    + (q2 .* exp(-k .* (n - 2))) .* (n >= 2);

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_series = 1:numel(combinedData)
    %Extraction of parameters
    q2 = mean(combinedData(i_series).separateQ2Vector);
    phi = mean(combinedData(i_series).separatePhiVector);
    k = mean(combinedData(i_series).separateKVector);

    %Calculation
    pspace = func_singlecellProb_MOIstar2(q2, phi, k, nSpace);

    %Visualization
    subplot(2, 3, i_series);
    hold on;

    plot(nSpace, pspace, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r',...
        'LineStyle', '-', 'Color', 'r');

    %Adjusting plot properties
    axis square;
    set(gca, 'YScale', 'log');
    xlim([0 10.5]);
    ylim([1e-4 1]);
    xlabel('Single-cell MOI');
    ylabel('Probability of lysogenization');
    title(['Growth rate = ' num2str(combinedData(i_series).growthRateInfected(1), '%.2f') ' per hour']);
    grid on; box on;
end

%% TRENDS OF FITTING PARAMETERS VS GROWTH RATE
clc
close all

%Some final preparations
for i_series = 1:numel(combinedData)
    combinedData(i_series).cellConcInfected = combinedData(i_series).odInfected * 1e8;
end

%Fitting functions
samplingN = 1000;
gSpace = 0:0.05:2;
linear.function = @(a, b, x) a.*x + b;
linear.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [-Inf -Inf], 'Upper', [Inf Inf], 'Startpoint', [0 0]);
linear.model = fittype(linear.function, 'Options', linear.options,...
    'Independent', 'x', 'Coefficients', {'a', 'b'});

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

growthRateVectorScatter = [];
growthRateVector = [];
for i_series = 1:numel(combinedData)
    growthRateVectorScatter = [growthRateVectorScatter,...
        combinedData(i_series).growthRateInfectedVector(1),...
        combinedData(i_series).growthRateInfectedVector(2)];

    growthRateVector = [growthRateVector combinedData(i_series).growthRateInfected(1)];
end

%LYSOGENIZATION PROBABILITY AT SC-MOI = 1
subplot(1, 2, 1);
hold on;

Q1VectorScatter = [];
Q1Vector = [];
Q1Vector_error = [];
for i_series = 1:numel(combinedData)
    Q1VectorScatter = [Q1VectorScatter combinedData(i_series).separateQ2Vector.*combinedData(i_series).separatePhiVector];

    Q1Vector = [Q1Vector mean(combinedData(i_series).separateQ2Vector.*combinedData(i_series).separatePhiVector)];
    Q1Vector_error = [Q1Vector_error std(combinedData(i_series).separateQ2Vector.*combinedData(i_series).separatePhiVector)/sqrt(numel(combinedData(i_series).separateQ2Vector))];
end

%With parametrization (Bootstrapped)
paramSampling = zeros(samplingN, 2);
for i_sampling = 1:samplingN
    [~, idx] = datasample(growthRateVectorScatter, numel(growthRateVectorScatter));

    [xToFit, yToFit] = prepareCurveData(growthRateVectorScatter(idx), log10(Q1VectorScatter(idx)));
    [param_q1, gof] = fit(xToFit, yToFit, linear.model);
    paramSampling(i_sampling, :) = [param_q1.a, param_q1.b];
end
paramMaster.q1 = mean(paramSampling, 1);
paramMaster.q1_SE = std(paramSampling, [], 1);
ySpace = linear.function(paramMaster.q1(1), paramMaster.q1(2), gSpace);
ySpace_propSE = sqrt((gSpace .* paramMaster.q1_SE(1)).^2 + (paramMaster.q1_SE(2)).^2);

%Visualization
patch([gSpace fliplr(gSpace)],...
    [(ySpace-ySpace_propSE) fliplr(ySpace+ySpace_propSE)],...
    'r', 'EdgeColor', 'None', 'FaceAlpha', 0.1);
graphParam = plot(gSpace, ySpace, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'r',...
    'DisplayName', ['Linear in log (Bootstrap SE):' newline,...
    '\beta_1 = ' num2str(paramMaster.q1(1), '%.2f') ' \pm ' num2str(paramMaster.q1_SE(1), '%.2f') newline,...
    '\beta_0 = ' num2str(paramMaster.q1(2), '%.2f') ' \pm ' num2str(paramMaster.q1_SE(2), '%.2f')]);
for i_repeat = 1:numel(runMarkers)
    graphScatter(i_repeat) = plot(growthRateVectorScatter(i_repeat:2:numel(growthRateVectorScatter)),...
        log10(Q1VectorScatter(i_repeat:2:numel(Q1VectorScatter))),...
        'Marker', runMarkers(i_repeat), 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k', 'LineStyle', 'None',...
        'DisplayName', ['Fitted values of Run #' num2str(i_repeat)]);
end
graphStats = errorbar(growthRateVector, log10(Q1Vector), sqrt((1./(Q1Vector .* log(10)) .* Q1Vector_error).^2),...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'LineStyle', 'None', 'Color', 'k',...
    'DisplayName', 'Mean \pm SE');

%Adjusting plot properties
axis square;
legend([graphParam], 'Location', 'SouthWest', 'FontSize', 10);
xlim([0 1.2]);
xlabel('Growth rate (/hr)');
ylim([-4.5 -0.5]);
ylabel('log_{10}(Q1)');
title(['Q1 vs. growth rate'], 'FontSize', 12);
grid on; box on;

%PHI, RATIO OF Q1 AND Q2
subplot(1, 2, 2);
hold on;

phiRecipVectorScatter = [];
phiRecipVector = [];
phiRecipVector_error = [];
for i_series = 1:numel(combinedData)
    phiRecipVectorScatter = [phiRecipVectorScatter 1./combinedData(i_series).separatePhiVector];

    phiRecipVector = [phiRecipVector mean(1./[combinedData(i_series).separatePhiVector])];
    phiRecipVector_error = [phiRecipVector_error std(1./[combinedData(i_series).separatePhiVector])/sqrt(numel(combinedData(i_series).separatePhiVector))];
end

%With parametrization
paramSampling = zeros(samplingN, 3);
phiRecipStats = mean(log10(phiRecipVectorScatter(growthRateVectorScatter >= 0.2)));
phiRecipStats_SE = std(log10(phiRecipVectorScatter(growthRateVectorScatter >= 0.2))) / sqrt(sum(growthRateVectorScatter >= 0.2));

%Visualization
patch([min(gSpace) max(gSpace) max(gSpace) min(gSpace)],...
    [(phiRecipStats-phiRecipStats_SE) (phiRecipStats-phiRecipStats_SE) fliplr(phiRecipStats+phiRecipStats_SE) fliplr(phiRecipStats+phiRecipStats_SE)],...
    'r', 'EdgeColor', 'None', 'FaceAlpha', 0.1);
graphParam = plot(gSpace, repelem(phiRecipStats, 2, numel(gSpace)),...
    'LineStyle', '-', 'LineWidth', 1.5, 'Color', 'r',...
    'DisplayName', ['Average value before stationary phase']);
for i_repeat = 1:numel(runMarkers)
    graphScatter(i_repeat) = plot(growthRateVectorScatter(i_repeat:2:numel(growthRateVectorScatter)),...
        log10(phiRecipVectorScatter(i_repeat:2:numel(phiRecipVectorScatter))),...
        'Marker', runMarkers(i_repeat), 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k', 'LineStyle', 'None',...
        'DisplayName', ['Fitted values of Run #' num2str(i_repeat)]);
end
graphStats = errorbar(growthRateVector, log10(phiRecipVector), sqrt((1./(phiRecipVector .* log(10)) .* phiRecipVector_error).^2),...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'LineStyle', 'None', 'Color', 'k',...
    'DisplayName', 'Mean \pm SE');

%Adjusting plot properties
axis square;
legend([graphParam(1)], 'Location', 'SouthEast', 'FontSize', 10);
xlim([0 1.2]);
xlabel('Growth rate (/hr)');
% ylim([0.5 1e3]);
ylabel('log_{10}Q2/Q1');
title(['Ratio of Q2 to Q1 vs. growth rate'], 'FontSize', 12);
grid on; box on;

sgtitle(['Simple guides-to-the-eye (Not parametrization for 2D surface)']);

%% PARAMETRIZATION FOR 2D SURFACE (No bootstrapping)
clc
close all

%Some final preparations
for i_series = 1:numel(combinedData)
    combinedData(i_series).cellConcInfected = combinedData(i_series).odInfected * 1e8;
end

%Fitting functions
samplingN = 100;
gSpace = 0:0.01:2;
linear.function = @(a, b, x) a.*x + b;
linear.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [-Inf -Inf], 'Upper', [Inf Inf], 'Startpoint', [0 0]);
linear.model = fittype(linear.function, 'Options', linear.options,...
    'Independent', 'x', 'Coefficients', {'a', 'b'});

piecewiseTwoSlopes.function = @(a, b, c, gstar, x) (a.*x + b) .* (x < gstar) + (c.*x + (a-c).*gstar + b) .* (x >= gstar);
piecewiseTwoSlopes.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [-Inf -Inf -Inf 0], 'Upper', [Inf Inf Inf 1.5], 'Startpoint', [0 0 0 0.9]);
piecewiseTwoSlopes.model = fittype(piecewiseTwoSlopes.function, 'Options', piecewiseTwoSlopes.options,...
    'Independent', 'x', 'Coefficients', {'a', 'b', 'c', 'gstar'});

piecewiseSlopeFirst.function = @(a, b, gstar, x) (a.*x + b) .* (x < gstar) + (a.*gstar + b) .* (x >= gstar);
piecewiseSlopeFirst.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [-Inf -Inf 0], 'Upper', [Inf Inf 2], 'Startpoint', [0 0 1.2]);
piecewiseSlopeFirst.model = fittype(piecewiseSlopeFirst.function, 'Options', piecewiseSlopeFirst.options,...
    'Independent', 'x', 'Coefficients', {'a', 'b', 'gstar'});

piecewiseFlatFirst.function = @(a, b, gstar, x) (a.*gstar + b) .* (x < gstar) + (a.*x + b) .* (x >= gstar);
piecewiseFlatFirst.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [-Inf -Inf 0], 'Upper', [Inf Inf 2], 'Startpoint', [0 0 1.2]);
piecewiseFlatFirst.model = fittype(piecewiseFlatFirst.function, 'Options', piecewiseFlatFirst.options,...
    'Independent', 'x', 'Coefficients', {'a', 'b', 'gstar'});

quadratic.function = @(a, b, c, x) a.*x.^2 + b.*x + c;
quadratic.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [-Inf -Inf -Inf], 'Upper', [Inf Inf Inf], 'Startpoint', [0 0 0]);
quadratic.model = fittype(quadratic.function, 'Options', quadratic.options,...
    'Independent', 'x', 'Coefficients', {'a', 'b', 'c'});

piecewiseQuadratic.function = @(g1, g2, x) (0) .* (x < g1)...
    - (x - g1) .* (x - g2) .* (x >= g1 & x < g2)...
    + (0) .* (x >= g2);
piecewiseQuadratic.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [0 1], 'Upper', [1 2], 'Startpoint', [0.7 1.6]);
piecewiseQuadratic.model = fittype(piecewiseQuadratic.function, 'Options', piecewiseQuadratic.options,...
    'Independent', 'x', 'Coefficients', {'g1', 'g2'});

%Preparations for Moldovan calculations
nBVector = logspace(6, 9, 50);
temp = 30; %Celsius
time = 30*60; %seconds of adsorption before dilution
k = 1.6e-11; %For 30C
ksingle = 1.9e-3; %For 30C
kdouble = 7.8e-4; %For 30C

tau1 = 1./(1/2 * ((k.*nBVector + ksingle + kdouble)...
    + sqrt((k.*nBVector + ksingle + kdouble).^2 - 4*k.*kdouble.*nBVector)));
tau2 = 1./(1/2 * ((k.*nBVector + ksingle + kdouble)...
    - sqrt((k.*nBVector + ksingle + kdouble).^2 - 4*k.*kdouble.*nBVector)));
fractionFree = 1 ./ (tau2 - tau1) .* ((1/kdouble - tau1) .* exp(-time./tau1)...
    - (1/kdouble - tau2) .* exp(-time./tau2));
fractionAdsorbed = 1 - fractionFree;

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

%X-SCALING PARAMETERS
subplot(2, 3, 1);
hold on;

cellConcVectorScatter = [];
growthRateVectorScatter = [];
xScaleVectorScatter = [];
cellConcVector = [];
growthRateVector = [];
xScaleVector = [];
xScaleVector_error = [];
for i_series = 1:numel(combinedData)
    xScaleVectorScatter = [xScaleVectorScatter combinedData(i_series).separateAVector];
    growthRateVectorScatter = [growthRateVectorScatter,...
        combinedData(i_series).growthRateInfectedVector(1),...
        combinedData(i_series).growthRateInfectedVector(2)];
    cellConcVectorScatter = [cellConcVectorScatter repelem(combinedData(i_series).cellConcInfected(1), 1, numel(combinedData(i_series).separateAVector))];

    xScaleVector = [xScaleVector mean(combinedData(i_series).separateAVector)];
    xScaleVector_error = [xScaleVector_error std(combinedData(i_series).separateAVector)/sqrt(numel(combinedData(i_series).separateAVector))];
    growthRateVector = [growthRateVector combinedData(i_series).growthRateInfected(1)];
    cellConcVector = [cellConcVector combinedData(i_series).cellConcInfected(1)];
end

%Visualization
graphMoldovan = plot(nBVector, fractionAdsorbed, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'r',...
    'DisplayName', ['Moldovan-predicted']);
graphScatter = plot(cellConcVectorScatter, xScaleVectorScatter,...
    'Marker', 'v', 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k', 'LineStyle', 'None',...
    'DisplayName', 'Fitted values');
graphStats = errorbar(cellConcVector, xScaleVector, xScaleVector_error,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'LineStyle', 'None', 'Color', 'k',...
    'DisplayName', 'Mean \pm SE');

%Adjusting plot properties
axis square;
legend([graphScatter graphMoldovan graphStats], 'Location', 'SouthEast', 'FontSize', 9);
xlim([0 6e8]);
xlabel('Cell concentration (CFU/mL)');
ylim([1e-2 5]);
ylabel('x-scaling factor');
title(['x-scaling factor vs. cell concentration'], 'FontSize', 12);
set(gca, 'YScale', 'Log');
grid on; box on;

subplot(2, 3, 2);
hold on;

%With parametrization (Bootstrapped)
paramSampling = zeros(samplingN, 3);
for i_sampling = 1:samplingN
    [~, idx] = datasample(growthRateVectorScatter, numel(growthRateVectorScatter));
%     idx = 1:numel(growthRateVectorScatter);

    [xToFit, yToFit] = prepareCurveData(growthRateVectorScatter(idx), log(xScaleVectorScatter(idx)));
    [param_xscale, gof] = polyfit(xToFit, yToFit, 2);
    paramSampling(i_sampling, :) = [param_xscale(1), param_xscale(2), param_xscale(3)];
end
paramMaster.xscale = mean(paramSampling, 1);
paramMaster.xscale_SE = std(paramSampling, [], 1);
ySpace = exp(paramMaster.xscale(1).*gSpace.^2 + paramMaster.xscale(2).*gSpace + paramMaster.xscale(3));

%Visualization
graphParam = plot(gSpace, ySpace, 'LineStyle', '-', 'LineWidth', 1, 'Color', 'r',...
    'DisplayName', ['Quadratic in log (Bootstrap SE):' newline,...
    '\beta_2 = ' num2str(paramMaster.xscale(1), '%.2f') ' \pm ' num2str(paramMaster.xscale_SE(1), '%.2f') newline,...
    '\beta_1 = ' num2str(paramMaster.xscale(2), '%.2f') ' \pm ' num2str(paramMaster.xscale_SE(2), '%.2f') newline,...
    '\beta_0 = ' num2str(paramMaster.xscale(3), '%.2f') ' \pm ' num2str(paramMaster.xscale_SE(3), '%.2f')]);
graphScatter = plot(growthRateVectorScatter, xScaleVectorScatter,...
    'Marker', 'v', 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k', 'LineStyle', 'None',...
    'DisplayName', 'Fitted values');
graphStats = errorbar(growthRateVector, xScaleVector, xScaleVector_error,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'LineStyle', 'None', 'Color', 'k',...
    'DisplayName', 'Mean \pm SE');

%Adjusting plot properties
axis square;
legend([graphParam], 'Location', 'SouthWest', 'FontSize', 9);
xlim([0 1.2]);
xlabel('Growth rate (/hr)');
ylim([1e-2 5]);
ylabel('x-scaling factor');
title(['vs. growth rate'], 'FontSize', 12);
set(gca, 'YScale', 'Log');
grid on; box on;

%REPRESSION FACTOR
subplot(2, 3, 3);
hold on;

kVectorScatter = [];
kVector = [];
kVector_error = [];
for i_series = 1:numel(combinedData)
    kVectorScatter = [kVectorScatter combinedData(i_series).separateKVector];

    kVector = [kVector mean(combinedData(i_series).separateKVector)];
    kVector_error = [kVector_error std(combinedData(i_series).separateKVector)/sqrt(numel(combinedData(i_series).separateKVector))];
end

%With parametrization (Bootstrapped)
paramSampling = zeros(samplingN, 2);
for i_sampling = 1:samplingN
    [~, idx] = datasample(growthRateVectorScatter, numel(growthRateVectorScatter));
%     idx = 1:numel(growthRateVectorScatter);

    [xToFit, yToFit] = prepareCurveData(growthRateVectorScatter(idx), kVectorScatter(idx));
    [param_k, gof] = fit(xToFit, yToFit, piecewiseQuadratic.model);
    paramSampling(i_sampling, :) = [param_k.g1, param_k.g2];
end
paramMaster.k = mean(paramSampling, 1);
paramMaster.k_SE = std(paramSampling, [], 1);
ySpace = piecewiseQuadratic.function(paramMaster.k(1), paramMaster.k(2), gSpace);

%Visualization
graphParam = plot(gSpace, ySpace, 'LineStyle', '-', 'LineWidth', 1, 'Color', 'r',...
    'DisplayName', ['Piecewise quadratic (Bootstrap SE):' newline,...
    'g_1 = ' num2str(paramMaster.k(1), '%.2f') ' \pm ' num2str(paramMaster.k_SE(1), '%.2f') newline,...
    'g_2 = ' num2str(paramMaster.k(2), '%.2f') ' \pm ' num2str(paramMaster.k_SE(2), '%.2f')]);
graphScatter = plot(growthRateVectorScatter, kVectorScatter,...
    'Marker', 'v', 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k', 'LineStyle', 'None',...
    'DisplayName', 'Fitted values');
graphStats = errorbar(growthRateVector, kVector, kVector_error,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'LineStyle', 'None', 'Color', 'k',...
    'DisplayName', 'Mean \pm SE');

%Adjusting plot properties
axis square;
legend([graphParam], 'Location', 'NorthWest', 'FontSize', 7);
xlim([0 1.2]);
xlabel('Growth rate (/hr)');
ylim([0 0.18]);
ylabel('Repression factor');
title(['Repression factor vs. growth rate'], 'FontSize', 12);
% set(gca, 'YScale', 'Log');
grid on; box on;

%LYSOGENIZATION PROBABILITY AT SC-MOI = 1
subplot(2, 3, 5);
hold on;

Q1VectorScatter = [];
Q1Vector = [];
Q1Vector_error = [];
for i_series = 1:numel(combinedData)
    Q1VectorScatter = [Q1VectorScatter combinedData(i_series).separateQ2Vector.*combinedData(i_series).separatePhiVector];

    Q1Vector = [Q1Vector mean(combinedData(i_series).separateQ2Vector.*combinedData(i_series).separatePhiVector)];
    Q1Vector_error = [Q1Vector_error std(combinedData(i_series).separateQ2Vector.*combinedData(i_series).separatePhiVector)/sqrt(numel(combinedData(i_series).separateQ2Vector))];
end

%With parametrization (Bootstrapped)
paramSampling = zeros(samplingN, 2);
for i_sampling = 1:samplingN
    [~, idx] = datasample(growthRateVectorScatter, numel(growthRateVectorScatter));
%     idx = 1:numel(growthRateVectorScatter);

    [xToFit, yToFit] = prepareCurveData(growthRateVectorScatter(idx), log(Q1VectorScatter(idx)));
    [param_q1, gof] = fit(xToFit, yToFit, linear.model);
    paramSampling(i_sampling, :) = [param_q1.a, param_q1.b];
end
paramMaster.q1 = mean(paramSampling, 1);
paramMaster.q1_SE = std(paramSampling, [], 1);
ySpace = exp(linear.function(paramMaster.q1(1), paramMaster.q1(2), gSpace));

%Visualization
graphParam = plot(gSpace, ySpace, 'LineStyle', '-', 'LineWidth', 1, 'Color', 'r',...
    'DisplayName', ['Linear in log (Bootstrap SE):' newline,...
    '\beta_1 = ' num2str(paramMaster.q1(1), '%.2f') ' \pm ' num2str(paramMaster.q1_SE(1), '%.2f') newline,...
    '\beta_0 = ' num2str(paramMaster.q1(2), '%.2f') ' \pm ' num2str(paramMaster.q1_SE(2), '%.2f')]);
graphScatter = plot(growthRateVectorScatter, Q1VectorScatter,...
    'Marker', 'v', 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k', 'LineStyle', 'None',...
    'DisplayName', 'Fitted values');
graphStats = errorbar(growthRateVector, Q1Vector, Q1Vector_error,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'LineStyle', 'None', 'Color', 'k',...
    'DisplayName', 'Mean \pm SE');

%Adjusting plot properties
axis square;
legend([graphParam], 'Location', 'SouthWest', 'FontSize', 9);
xlim([0 1.2]);
xlabel('Growth rate (/hr)');
ylim([1e-5 1]);
ylabel('Q1');
title(['Q1 vs. growth rate'], 'FontSize', 12);
set(gca, 'YScale', 'Log');
grid on; box on;

%LYSOGENIZATION PROBABILITY AT SC-MOI = 2
subplot(2, 3, 6);
hold on;

Q2VectorScatter = [];
Q2Vector = [];
Q2Vector_error = [];
for i_series = 1:numel(combinedData)
    Q2VectorScatter = [Q2VectorScatter combinedData(i_series).separateQ2Vector];

    Q2Vector = [Q2Vector mean(combinedData(i_series).separateQ2Vector)];
    Q2Vector_error = [Q2Vector_error std(combinedData(i_series).separateQ2Vector)/sqrt(numel(combinedData(i_series).separateQ2Vector))];
end

%With parametrization (Bootstrapped)
paramSampling = zeros(samplingN, 4);
for i_sampling = 1:samplingN
    [~, idx] = datasample(growthRateVectorScatter, numel(growthRateVectorScatter));
    % idx = 1:numel(growthRateVectorScatter);

    [xToFit, yToFit] = prepareCurveData(growthRateVectorScatter(idx), log(Q2VectorScatter(idx)));
    [param_q2, gof] = fit(xToFit, yToFit, piecewiseTwoSlopes.model);
    paramSampling(i_sampling, :) = [param_q2.a, param_q2.b, param_q2.c, param_q2.gstar];
end
paramMaster.q2 = mean(paramSampling, 1);
paramMaster.q2_SE = std(paramSampling, [], 1);
ySpace = exp(piecewiseTwoSlopes.function(paramMaster.q2(1), paramMaster.q2(2), paramMaster.q2(3), paramMaster.q2(4), gSpace));

%Visualization
graphParam = plot(gSpace, ySpace, 'LineStyle', '-', 'LineWidth', 1, 'Color', 'r',...
    'DisplayName', ['Piecewise linear in log (Bootstrap SE):' newline,...
    '\beta_1 = ' num2str(paramMaster.q2(1), '%.2f') ' \pm ' num2str(paramMaster.q2_SE(1), '%.2f') newline,...
    '\beta_0 = ' num2str(paramMaster.q2(2), '%.2f') ' \pm ' num2str(paramMaster.q2_SE(2), '%.2f') newline,...
    '\beta_2 = ' num2str(paramMaster.q2(3), '%.2f') ' \pm ' num2str(paramMaster.q2_SE(2), '%.2f') newline,...
    'g* = ' num2str(paramMaster.q2(4), '%.2f') ' \pm ' num2str(paramMaster.q2_SE(3), '%.2f')]);
graphScatter = plot(growthRateVectorScatter, Q2VectorScatter,...
    'Marker', 'v', 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k', 'LineStyle', 'None',...
    'DisplayName', 'Fitted values');
graphStats = errorbar(growthRateVector, Q2Vector, Q2Vector_error,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None', 'LineStyle', 'None', 'Color', 'k',...
    'DisplayName', 'Mean \pm SE');

%Adjusting plot properties
axis square;
legend([graphParam], 'Location', 'SouthWest', 'FontSize', 9);
xlim([0 1.2]);
xlabel('Growth rate (/hr)');
ylim([1e-5 1]);
ylabel('Q2');
title(['Q2 vs. growth rate'], 'FontSize', 12);
set(gca, 'YScale', 'Log');
grid on; box on;

sgtitle('Ad hoc parametrization for 2D surface later');

%% 2D RESPONSE SURFACE - DATA INTERPOLATION VS. MODEL PREDICTIONS
clc
close all

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.15 0.9 0.7]);
hold on;

%Custom color map
colormap parula;
xRange = [3e-3 50];
yRange = [0 1.25];
zRange = [-4 0];

%INTERPOLATION OF RAW DATA
subplot(1, 3, 1);
hold on;

%Preparations
nPoints = 30;
moifor3D = [];
growthRatefor3D = [];
fractionfor3D = [];

for i_series = 1:numel(data)
    %Data extraction
    i_run = data(i_series).run;

    %Correcting for the x-scaling factor
    currentG = combinedData(mod(i_series-1, 6)+1).growthRateInfected(1);
    xScale = exp(quadratic.function(paramMaster.xscale(1), paramMaster.xscale(2), paramMaster.xscale(3), currentG));

    MOIseries = [data(i_series).MOIseries];
    growthRate = [data(i_series).growthRateInfected];
    fraction = [data(i_series).fractionLysogeny];
    flagValid = fraction(1, :) > 0;

    moifor3D = [moifor3D MOIseries(flagValid).*xScale];
    growthRatefor3D = [growthRatefor3D repelem(growthRate(1, :), 1, sum(flagValid))];
    fractionfor3D = [fractionfor3D log10(fraction(1, flagValid))];
end

%Extrapolation
moiV = logspace(log10(min(moifor3D)), log10(max(moifor3D)), nPoints);
growthRateV = linspace(min(growthRatefor3D), max(growthRatefor3D), nPoints);
[moiMesh, growthRateMesh] = meshgrid(moiV, growthRateV);
fractionMesh = griddata(moifor3D, growthRatefor3D, fractionfor3D, moiMesh, growthRateMesh,...
    'natural');

%Visualization of the interpolation surface
graphSurf = surf(moiMesh, growthRateMesh, fractionMesh,...
    'EdgeColor', 'None',...
    'DisplayName', 'Interpolation');

%And of individual replicate runs
for i_series = 1:numel(data)
    %Data extraction
    i_run = data(i_series).run;

    %Correcting for the x-scaling factor
    currentG = combinedData(mod(i_series-1, 6)+1).growthRateInfected(1);
    xScale = exp(quadratic.function(paramMaster.xscale(1), paramMaster.xscale(2), paramMaster.xscale(3), currentG));

    MOIseries = [data(i_series).MOIseries];
    growthRate = [data(i_series).growthRateInfected];
    fraction = [data(i_series).fractionLysogeny];
    flagValid = fraction(1, :) > 0;

    %Visualization
    graphData(i_series) = scatter3(MOIseries(flagValid).*xScale(1),...
        repelem(growthRate(1, :), 1, sum(flagValid)),...
        log10(fraction(1, flagValid)),...
        40, 'k', 'filled', 'Marker', runMarkers(i_run), 'MarkerEdgeColor', 'None',...
        'DisplayName', ['Run #' num2str(ceil(i_series/6))]);
end

%Adjusting plot properties
colormap(parula);
cbar = colorbar();
caxis(zRange);
set(cbar, 'Location', 'SouthOutside', 'Xtick', [-4 -3 -2 -1 0]);
cbar.Label.String = 'log_{10}(Fraction of lysogeny)';
cbar.Label.FontSize = 12;
% set(graphSurf, 'FaceAlpha', 0.8, 'FaceLighting', 'gouraud');
% view(200, 35);

axis square;
% legend([graphData([1 7]) graphSurf], 'Location', 'NorthWest', 'FontSize', 8);
xlim([1e-3 100]);
xlabel('Bulk MOI (Corrected)');
ylim(yRange);
ylabel('Growth rate (db/hr)');
% zlim([-6 0]);
zlabel('log10(Fraction of lysogeny)');
title('Interpolation from data', 'FontSize', 14);
set(gca, 'XScale', 'Log');
grid on; box on;

%WITH MEDIAN FILTERING IN 2D
subplot(1, 3, 2);
hold on;

%Preparations
moifor3D = [];
growthRatefor3D = [];
fractionfor3D = [];

for i_series = 1:numel(data)
    %Data extraction
    i_run = data(i_series).run;

    %Correcting for the x-scaling factor
    currentG = combinedData(mod(i_series-1, 6)+1).growthRateInfected(1);
    xScale = exp(quadratic.function(paramMaster.xscale(1), paramMaster.xscale(2), paramMaster.xscale(3), currentG));

    MOIseries = [data(i_series).MOIseries];
    growthRate = [data(i_series).growthRateInfected];
    fraction = [data(i_series).fractionLysogeny];
    flagValid = fraction(1, :) > 0;

    moifor3D = [moifor3D MOIseries(flagValid).*xScale];
    growthRatefor3D = [growthRatefor3D repelem(growthRate(1, :), 1, sum(flagValid))];
    fractionfor3D = [fractionfor3D log10(fraction(1, flagValid))];
end

%Extrapolation with filtering
moiV = logspace(log10(min(moifor3D)), log10(max(moifor3D)), nPoints);
growthRateV = linspace(min(growthRatefor3D), max(growthRatefor3D), nPoints);
[moiMesh, growthRateMesh] = meshgrid(moiV, growthRateV);
fractionMesh = griddata(moifor3D, growthRatefor3D, fractionfor3D, moiMesh, growthRateMesh,...
    'natural');
filteredFractionMesh = medfilt2(fractionMesh, [2 5], 'symmetric');

%Visualization of the interpolation surface
graphSurf = surf(moiMesh, growthRateMesh, filteredFractionMesh,...
    'EdgeColor', 'None',...
    'DisplayName', 'Interpolation');

%Adjusting plot properties
colormap(parula);
cbar = colorbar();
caxis(zRange);
set(cbar, 'Location', 'SouthOutside', 'Xtick', [-4 -3 -2 -1 0]);
cbar.Label.String = 'log_{10}(Fraction of lysogeny)';
cbar.Label.FontSize = 12;
% set(graphSurf, 'FaceAlpha', 0.8, 'FaceLighting', 'gouraud');
% view(200, 35);

axis square;
% legend([graphData([1 7]) graphSurf], 'Location', 'NorthWest', 'FontSize', 8);
xlim(xRange);
xlabel('Bulk MOI (Corrected)');
ylim(yRange);
ylabel('Growth rate (db/hr)');
xticks([0.001 0.01 0.1 1 10]);
yticks([0:0.25:1]);
% zlim([-6 0]);
zlabel('log10(Fraction of lysogeny)');
title(['With 2D median filtering'], 'FontSize', 14);
set(gca, 'XScale', 'Log');
grid on; box on;

%PARAMETRIZED SURFACE
subplot(1, 3, 3);
hold on;

%Surface
moiV = logspace(log10(xRange(1)), log10(xRange(2)), nPoints);
growthRateV = linspace(yRange(1), yRange(2), nPoints);
[moiMesh, growthRateMesh] = meshgrid(moiV, growthRateV);
fractionMesh = size(moiMesh);
params = [paramMaster.xscale paramMaster.k paramMaster.q1 paramMaster.q2];
for i_moi = 1:numel(moiV)
    currentM = moiV(i_moi);
    for i_g = 1:numel(growthRateV)
        currentG = growthRateV(i_g);

        %Calculations
        currentFraction = exp(func_gParametrizedLysogeny_noshift(params, [currentG, currentM]));
        fractionMesh(i_g, i_moi) = log10(currentFraction);
    end
end
 
%Visualization
graphSurf = surf(moiMesh, growthRateMesh, fractionMesh,...
    'EdgeColor', 'None',...
    'DisplayName', 'Predicted surface');

%Adjusting plot properties
colormap(parula);
cbar = colorbar();
caxis(zRange);
set(cbar, 'Location', 'SouthOutside', 'Xtick', [-4 -3 -2 -1 0]);
cbar.Label.String = 'log_{10}(Fraction of lysogeny)';
cbar.Label.FontSize = 12;
% set(graphSurf, 'FaceAlpha', 0.8, 'FaceLighting', 'gouraud');
% view(200, 35);

axis square;
% legend([graphData([1 7]) graphSurf], 'Location', 'NorthWest', 'FontSize', 8);
xlim(xRange);
xlabel('Bulk MOI');
ylim(yRange);
ylabel('Growth rate (db/hr)');
xticks([0.001 0.01 0.1 1 10]);
yticks([0:0.25:1]);
% zlim([-6 0]);
zlabel('log10(Fraction of lysogeny)');
title('Model predictions', 'FontSize', 14);
set(gca, 'XScale', 'Log');
grid on; box on;

sgtitle('Predictions using ad hoc parametrization of Q1 and Q2');