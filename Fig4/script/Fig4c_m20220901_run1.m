close all
clear all
clc

%INFECTION PLATE
%Data loading and processing
folderPath = [''];
fileName = ['20220901_ODbasedMOIresponse_optimalphysio_infectionplate.xlsx'];

%Data processing
timeVectorInfection = readmatrix([folderPath fileName], 'Range', 'B37:J37')/60;
rawInfectionOD = readmatrix([folderPath fileName], 'Range', 'B39:J86');

%Filling in the infection replicates -- See the plate map for more details
currentCopy = rawInfectionOD(1, :);
for i_well = 1:size(rawInfectionOD, 1)
    if i_well <= 42
        replicateTmp = 3;
    else
        replicateTmp = 2;
    end
    if mod(i_well, replicateTmp) == 1
        currentCopy = rawInfectionOD(i_well, :);
    else
        rawInfectionOD(i_well, :) = currentCopy;
    end
end

idBlank = [47 48];
blankOD = mean(rawInfectionOD(idBlank, :), 1);
blankedInfectionOD = rawInfectionOD - blankOD;

%DETECTION PLATE
%Data loading and processing
fileName = ['20220901_ODbasedMOIresponse_optimalphysio_detectionplate.xlsx'];

%Data processing
timeVectorDetection = readmatrix([folderPath fileName], 'Range', 'B37:GQ37')/60;
rawDetectionOD = readmatrix([folderPath fileName], 'Range', 'B39:GQ86');

idBlank = [48];
blankOD = mean(rawDetectionOD(idBlank, :), 1);
blankedDetectionOD = rawDetectionOD - blankOD;

%Concatenating the two plates together
timeVector = [timeVectorInfection - (timeVectorInfection(end) + 1), timeVectorDetection];
rawOD = [rawInfectionOD, rawDetectionOD];
blankedOD = [blankedInfectionOD, blankedDetectionOD];

disp('Data import and processing complete.');

%Biological information of the samples
conditionName = ["Infection of \lambda P+ in LBM at 30^oC, with optimal-physiology MG1655. Run 1."];

cellConcentration = [4.68e7]; %CFU/mL
cellVolume = 500; %uL
phageStockConcentration = 2.5e11; %PFU/mL
phageVolume = 10; %uL
MOIfolds = 2;
MOIseries = (phageStockConcentration .* MOIfolds.^(0:-1:-13) .* phageVolume) ./ (cellConcentration .* cellVolume);
% oneCellOD = 1/(2e9);

%Technical information for the TECAN plate
sampleReplicates = [repelem(3, 1, 14) 2 2 2];
nSamples = numel(sampleReplicates);
sampleNames = cell(1, nSamples); %The assignment method below only works for cells
for i_sample = 1:nSamples
    if ismember(i_sample, [1:14])
        sampleNames{i_sample} = strcat("MOI = ", string(num2str(MOIseries(i_sample), '%.2e')));
    elseif i_sample == 15
        sampleNames{i_sample} = "Uninfected (+Kan)";
    elseif i_sample == 16
        sampleNames{i_sample} = "Uninfected (-Kan)";
    elseif i_sample == 17
        sampleNames{i_sample} = "Blank LBM";
    end
end
sampleNames = string(sampleNames); %Converting from cell to array for easy access
sampleColors = [hsv(size(MOIseries, 2))*0.7; repmat([1 1 1]*0.3, 3, 1)];
sampleLines = [repelem("-", 1, size(MOIseries, 2)) "--" "-" ":"];

wellMap = repelem(1:nSamples, 1, sampleReplicates); %Which well is which sample
nWells = numel(wellMap);
[~, firstWell, ~] = unique(wellMap);
firstWell = firstWell';
wellNames = repelem(sampleNames, 1, sampleReplicates);
wellColors = repelem(sampleColors, sampleReplicates, 1);
wellLines = repelem(sampleLines, 1, sampleReplicates);

%Further instructions for plotting
nCol = 4;
nRow = 4;
xLimit = 1200; %minutes

disp('Instructions for plotting complete.');

%% OD600 curves of all samples
close all
clc

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_type = 1:2
    subplot(1, 2, i_type);
    hold on;

    switch i_type
        case 1
            currentOD = rawOD;

            title(['Raw OD_{600} of both the Infection and Detection plates'], 'FontSize', 14);
            ylabel('OD_{600}');
            ylim([1e-4 2]);
        case 2
            currentOD = blankedOD;

            title(['Blank-subtracted OD_{600} of both plates'], 'FontSize', 14);
            ylabel('Blank-subtracted OD_{600}');
            ylim([1e-4 2]);
    end

    for i_well = 1:nWells
        graph(i_well) = plot(timeVector, currentOD(i_well, :),...
            'LineStyle', wellLines(i_well), 'LineWidth', 1.5,...
            'Color', wellColors(i_well, :),...
            'DisplayName', wellNames(i_well));
    end

    xline(timeVector(6), '--k');
    xline(timeVector(9), '--k');
    xline(timeVector(15), '--k');

    %Adjusting plot properties
    legend([graph(firstWell)],...
        'Location', 'SouthEast', 'FontSize', 10);
    xlim([-60 xLimit]);
    xlabel('Time (min)');
    set(gca, 'YScale', 'Log');
    grid on; box on;
end

sgtitle([conditionName]);

%% Fitting for growth rates
clc
close all

figure('Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.9 0.9]);
hold on;

%Final preparations for fitting
minTime = 5; %minutes
rangeOD = [0.01 0.1]; %in blank-subtracted OD values, between which to fit
xspace = 0:1:(xLimit);

generationTimeCell = cell(1, nSamples);
generationTimeArray = zeros(3, nSamples);

initialOD = zeros(3, nSamples);
initialODCell = cell(1, nSamples);

for i_sample = 1:(nSamples-1)
    subplot(nRow, nCol, i_sample);
    hold on;

    %Range for growth rate fitting
    yline(rangeOD(1), '--k');
    yline(rangeOD(2), '--k');

    relevantWells = find((wellMap == i_sample));

    gVector = [];
    initialVector = [];
    for i_well = relevantWells
        %Fitting for growth rates
        flagRange = (timeVector >= minTime) & (blankedOD(i_well, :) >= rangeOD(1)) & (blankedOD(i_well, :) <= rangeOD(2));

        if sum(flagRange) >= 3
            param = polyfit(timeVector(flagRange), log(blankedOD(i_well, flagRange)), 1);
            yspace = exp(param(1) .* xspace + param(2));

            %Saving for later
            gVector = [gVector log(2)/param(1)];
            initialVector = [initialVector exp(param(2))];

            graphFit(i_well) = plot(xspace, yspace,...
                'LineStyle', '-', 'Color', 'k',...
                'DisplayName', ['Fitting for ' char(wellNames(i_well))]);
        else
            gVector = [gVector NaN]; %NaN for no detectable growth
            initialVector = [initialVector 0]; %Detection limit of one cell, if no growth
        end

        %Visualization
        graphData(i_well) = plot(timeVector, blankedOD(i_well, :),...
            'LineStyle', ':', 'LineWidth', 1.5,...
            'Color', wellColors(i_well, :),...
            'DisplayName', wellNames(i_well));
    end

    %Statistics calculation
    generationTimeCell{i_sample} = gVector;

    generationTimeArray(1, i_sample) = mean(gVector, 'omitnan');
    generationTimeArray(2, i_sample) = std(gVector, 'omitnan');
    generationTimeArray(3, i_sample) = std(gVector)/sqrt(numel(gVector));

    initialODCell{i_sample} = initialVector;

    initialOD(1, i_sample) = mean(initialVector, 'omitnan');
    initialOD(2, i_sample) = std(initialVector, 'omitnan');
    initialOD(3, i_sample) = std(initialVector)/sqrt(numel(initialVector));

    %Adjusting plot properties
    xlim([0 xLimit]);
    xlabel('Time (min)');
    ylim([1e-4 3]);
    ylabel('OD_{600}');
    title(sampleNames(i_sample), 'FontSize', 12);
    set(gca, 'YScale', 'Log');
    grid on; box on;
end

sgtitle([conditionName]);

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_sample = 1:(nSamples-1)
    bar(i_sample, generationTimeArray(1, i_sample),...
        'FaceColor', sampleColors(i_sample, :));

    errorbar(i_sample, generationTimeArray(1, i_sample), generationTimeArray(3, i_sample),...
        'Marker', 'None', 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1.2, 'CapSize', 10);
end

%Adjusting plot properties
ylim([0 60]);
ylabel('Generation time (min)');
xticks(1:(nSamples-1));
xticklabels(sampleNames(1:(nSamples-1)));
title([char(conditionName) ' Generation time of all samples'], 'FontSize', 15);
grid on; box on

%% By extrapolation to the initial OD, for the relative fraction of lysogens
clc
close all

%Meta information
referenceSample = 16; %Uninfected without Kan
individualFractions = cell(1, nSamples);
relativeFraction = zeros(2, numel(MOIseries));

%Calculation of the lysogenization frequency
for i_sample = 1:numel(MOIseries)
    %From readings of individual replicates; the uninfected sample is still using statistics
    individualFractions{i_sample} = initialODCell{i_sample} / initialOD(1, referenceSample);
    relativeFraction(1, i_sample) = mean(individualFractions{i_sample});
    relativeFraction(2, i_sample) = std(individualFractions{i_sample})/sqrt(numel(individualFractions{i_sample}));

    %When mean >= std, for display purpose
    if relativeFraction(2, i_sample) >= relativeFraction(1, i_sample)
        relativeFraction(2, i_sample) = relativeFraction(2, i_sample) * (1 - 1e-4);
    end
end

moiReplicates = [];
fractionReplicates = [];
for i_moi = 1:numel(MOIseries)
    moiReplicates = [moiReplicates repelem(MOIseries(i_moi), 1, numel(individualFractions{i_moi}))];
    fractionReplicates = [fractionReplicates individualFractions{i_moi}];
end

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);

subplot(1, 2, 1);
hold on;

for i_sample = 1:numel(initialODCell)
    scatter(repelem(i_sample, 1, numel(initialODCell{i_sample})), initialODCell{i_sample},...
        'Marker', 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k',...
        'DisplayName', ['Of individual detection replicates']);
end

errorbar(1:numel(initialODCell), initialOD(1, :), initialOD(3, :),...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
    'LineStyle', 'None', 'Color', 'k', 'LineWidth', 1.2, 'CapSize', 8);

%Adjusting plot properties
xlim([0 numel(initialODCell)+1]);
xticks(1:numel(initialODCell));
xticklabels(sampleNames);
xlabel('Samples');
ylabel('Initial cell density (OD unit)');
title(['Estimation of the initial cell density using the extrapolation method'], 'FontSize', 14);
set(gca, 'YScale', 'Log');
grid on; box on;

subplot(1, 2, 2);
hold on;
axis square;

graphReplicates = scatter(moiReplicates, fractionReplicates,...
    'Marker', 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k',...
    'DisplayName', ['Of individual detection replicates']);

graphOD = errorbar(MOIseries(:), relativeFraction(1, :), relativeFraction(2, :),...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
    'LineStyle', 'None', 'Color', 'k', 'LineWidth', 1.2, 'CapSize', 8,...
    'DisplayName', ['By optical density, under Kan50 selection']);

%Adjusting plot properties
axis([1e-2 2e2 1e-6 1e-0]);
legend([graphOD graphReplicates], 'Location', 'NorthWest', 'FontSize', 11);
set(gca, 'XScale', 'Log', 'YScale', 'Log');
xlabel('Mixing MOI');
ylabel('Fraction of lysogeny');
title(['Fraction of lysogeny, without fitting'], 'FontSize', 14);
grid on; box on;

sgtitle(conditionName);

%% Modeling: Weitz multiple-step-functions, fitting in logarithmic space
clc
close all

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

%Defining the fit models
xspace = logspace(-5, 4, 100);
weitz(1).options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [0 0],...
    'Upper', [1 1],...
    'Startpoint', [0.1 1]);
weitz(1).function = @(a, q1, M)...
    log(q1) + log(1 - exp(-a.*M));
weitz(1).model = fittype(weitz(1).function,...
    'Options', weitz(1).options, 'Independent', 'M', 'Coefficient', {'a', 'q1'});

weitz(2).options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [0 0 0],...
    'Upper', [1 1 1],...
    'Startpoint', [0.1 0 1]);
weitz(2).function = @(a, q1, q2, M)...
    log(a.*M.*exp(-a.*M) .* q1 + (1 - exp(-a.*M) - a.*M.*exp(-a.*M)) .* (q1 + q2));
weitz(2).model = fittype(weitz(2).function,...
    'Options', weitz(2).options, 'Independent', 'M', 'Coefficient', {'a', 'q1', 'q2'});

weitz(3).options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [0 0 0 0],...
    'Upper', [1 1 1 1],...
    'Startpoint', [0.1 1 1 1]);
weitz(3).function = @(a, q1, q2, q3, M)...
    log(a.*M.*exp(-a.*M) .* q1 + (a.*M).^2./2.*exp(-a.*M) .* (q1 + q2) + (1 - exp(-a.*M) - a.*M.*exp(-a.*M) - (a.*M).^2./2.*exp(-a.*M)) .* (q1 + q2 + q3));
weitz(3).model = fittype(weitz(3).function,...
    'Options', weitz(3).options, 'Independent', 'M', 'Coefficient', {'a', 'q1', 'q2', 'q3'});

%Meta instructions for plotting
maxN = 5;
flagValid = ~isnan(fractionReplicates);

for i_model = 1:1:numel(weitz)
    %FITTING TO POP-AVERAGED LYSOGENIZATION FREQUENCY
    subplot(2, numel(weitz), i_model);
    hold on;

    %Fitting
    [xDataToFit, yDataToFit] = prepareCurveData(moiReplicates(flagValid), log(fractionReplicates(flagValid)));
    [param, gof] = fit(xDataToFit, yDataToFit, weitz(i_model).model);
    confidence = confint(param);
    switch i_model
        case 1
            yspace = weitz(i_model).function(param.a, param.q1, xspace);
            maxQ = param.q1;
        case 2
            yspace = weitz(i_model).function(param.a, param.q1, param.q2, xspace);
            maxQ = param.q1 + param.q2;
        case 3
            yspace = weitz(i_model).function(param.a, param.q1, param.q2, param.q3, xspace);
            maxQ = param.q1 + param.q2 + param.q3;
    end

    %Visualization
    graphReplicates = scatter(moiReplicates, fractionReplicates,...
        18, 'Marker', 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k',...
        'DisplayName', ['Of individual reading replicates']);

    graphModel = plot(xspace, exp(yspace),...
        'LineStyle', '-', 'LineWidth', 2, 'Color', 'r',...
        'DisplayName', ['x-scale (a) = ' num2str(param.a, '%.2f') newline,...
        'RMSE = ' num2str(gof.rmse, '%.3f')]);

    graphData = errorbar(MOIseries, relativeFraction(1, :), relativeFraction(2, :),...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
        'LineStyle', 'None', 'Color', 'k', 'LineWidth', 1, 'CapSize', 8,...
        'DisplayName', ['Data by plating']);

    %Adjusting plot properties
    axis([1e-3 5e2 1e-6 2e0]);
    legend([graphData graphModel], 'Location', 'NorthWest', 'FontSize', 10);
    set(gca, 'XScale', 'Log', 'YScale', 'Log');
    xlabel('Bulk MOI (Actual)');
    ylabel('Fraction of lysogeny, f');
    title(['Multiple-step-function, fitted to MOI* = ' num2str(i_model)], 'FontSize', 14);
    grid on; box on;

    %SINGLE-CELL LYSOGENIZATION PROBABILITY
    subplot(2, numel(weitz), i_model+numel(weitz));
    hold on;

    %Preparations
    nspace = 0:1:maxN;
    switch i_model
        case 1
            qspace = [0 repelem(param.q1, 1, maxN-i_model+1)];
        case 2
            qspace = [0 param.q1 repelem(param.q1 + param.q2, 1, maxN-i_model+1)];
        case 3
            qspace = [0 param.q1 (param.q1 + param.q2) repelem(param.q1 + param.q2 + param.q3, 1, maxN-i_model+1)];
    end

    %Visualization
    stairs(nspace, qspace,...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
        'LineStyle', '-', 'LineWidth', 2.5, 'Color', [1 1 1]*0.6);

    %Adjusting plot properties
    legend(['Q_{max} = ' num2str(max(qspace), '%.4f')],...
        'Location', 'NorthEast', 'FontSize', 11);
    axis([0 maxN 0 0.5]);
    xlabel('Single-cell MOI, n');
    ylabel('Probability, Q_n');
    title(['Single-cell probability of lysogenization'], 'FontSize', 14);
    grid on; box on;
    %     set(gca, 'YScale', 'Log');
end

sgtitle([char(conditionName) ' Fitting is performed in logarithmic space.']);

%% Other visualizations of the model fitting, with MOI shifted, with probability of lysogenization normalized
clc
close all

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

%Preparations
i_model = 2; %Best-fit variant
maxN = 5;
xspace = logspace(-4, 4, 100);
flagValidStats = ~isnan(relativeFraction(1, :));
flagValid = ~isnan(fractionReplicates);

%Fitting
[xDataToFit, yDataToFit] = prepareCurveData(moiReplicates(flagValid), log(fractionReplicates(flagValid)));
[param, gof] = fit(xDataToFit, yDataToFit, weitz(i_model).model);

switch i_model
    case 1
        yspace = weitz(i_model).function(param.a, param.q1, xspace);
        maxQ = param.q1;
    case 2
        yspace = weitz(i_model).function(param.a, param.q1, param.q2, xspace);
        maxQ = param.q1 + param.q2;
    case 3
        yspace = weitz(i_model).function(param.a, param.q1, param.q2, param.q3, xspace);
        maxQ = param.q1 + param.q2 + param.q3;
end

%Fraction of lysogeny vs. Average MOI
for i_type = 1:3
    subplot(2, 3, i_type);
    hold on;

    switch i_type
        case 1
            xScale = 1;
            yScale = 1;

            xlabel('Bulk MOI (Actual)');
            ylabel('Fraction of lysogeny, f');
            title(['Multiple-step-function, fitted to MOI* = ' num2str(i_model)], 'FontSize', 12);
        case 2
            xScale = param.a;
            yScale = 1;

            xlabel('Bulk MOI (Shifted)');
            ylabel('Fraction of lysogeny, f');
        case 3
            xScale = param.a;
            yScale = maxQ;

            xlabel('Bulk MOI (Shifted)');
            ylabel('Normalized fraction of lysogeny, f');
    end

    %Visualization
    graphModel = plot(xspace*xScale, exp(yspace)/yScale,...
        'LineStyle', '-', 'LineWidth', 2, 'Color', 'r',...
        'DisplayName', ['x-scale (a) = ' num2str(param.a, '%.2f') newline,...
        'RMSE = ' num2str(gof.rmse, '%.3f')]);

    graphData = errorbar(MOIseries(flagValidStats)*xScale, relativeFraction(1, flagValidStats)/yScale,...
        relativeFraction(2, flagValidStats)/yScale,...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
        'LineStyle', 'None', 'Color', 'k', 'LineWidth', 1, 'CapSize', 8);

    %Adjusting plot properties
    if i_type == 1
        legend(graphModel, 'Location', 'NorthWest', 'FontSize', 11);
    end
    axis([1e-4 5e2 1e-6 5e0]);
    set(gca, 'XScale', 'Log', 'YScale', 'Log');
    grid on; box on;
end

%Single-cell probability of lysogenization
%Preparations
nspace = 0:1:maxN;
switch i_model
    case 1
        qspace = [0 repelem(param.q1, 1, maxN-i_model+1)];
    case 2
        qspace = [0 param.q1 repelem(param.q1 + param.q2, 1, maxN-i_model+1)];
    case 3
        qspace = [0 param.q1 (param.q1 + param.q2) repelem(param.q1 + param.q2 + param.q3, 1, maxN-i_model+1)];
end

for i_type = [4 6]
    subplot(2, 3, i_type);
    hold on;
    
    switch i_type
        case 4
            yScale = 1;

            axis([0 maxN 0 maxQ*1.1]);
            xlabel('Single-cell MOI, n');
            ylabel('Probability, Q_n');
            title(['Single-cell probability of lysogenization'], 'FontSize', 12);
        case 6
            yScale = maxQ;

            axis([0 maxN 0 1.1]);
            xlabel('Single-cell MOI, n');
            ylabel('Normalized probability, Q_n');
    end

    %Visualization
    stairs(nspace, qspace/yScale,...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
        'LineStyle', '-', 'LineWidth', 2.5, 'Color', [1 1 1]*0.6);

    %Adjusting plot properties
    if i_type == 4
        legend(['Q_{max} = ' num2str(max(qspace), '%.4f')],...
            'Location', 'NorthWest', 'FontSize', 11);
    end
    grid on; box on;
end

sgtitle([conditionName]);