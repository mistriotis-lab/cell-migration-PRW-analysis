function [] = single_file_msd_process()
    %% single_file_msd_process
    % Takes a single cell trajectory data file for post processing.
    % Computes speed and persistence time from the mean squared
    % displacement (MSD) data fitted in a persistent random walk (PRW) 
    % model. Additionally, fits the exponent and anamolous diffusivity 
    % coefficient from the power law definition of MSD. 
    %
    % Lines 301 through 303 may be modified to change boundaries and
    % starting points of solved parameters in the fitted model.
    %
    %   1. begin by checking if track length is less than min_length time 
    %   steps. If so, do not analyze the track (insufficient data)
    %   
    %   2. Calculate tau values and the MSD of the series for that cell
    %
    %   3. Fit PRW model usings taus and calculated MSD values [using
    %   provided code from Wu et al., see function signature for full
    %   citation]
    %
    %   4. Goodness of Fit using R^2 and a fit flag for quick visual
    %   parsing
    %

clear all

tic

% Minimum track length required for analysis
min_length = 6;
% Target R^2 value to assign a status of "Good" fit
target_r2_val = 0.88;

% Get path to the data file being analyzed
[table_file_name, table_path] = uigetfile('*');
full_table_path = [table_path, table_file_name];
data = readtable(full_table_path);

% Check if time column label matches anticipated format
% If not, the label is changed within the data structure
if ismember('t_mins_', data.Properties.VariableNames)
        data.Properties.VariableNames{'t_mins_'} = 't_min_';
end

%% Pathing
% Path to an output directory within the data directory
outDir = fullfile(table_path, 'temp_out');
% Check if outDir exists, create folder if it does not
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

% Use fileparts to strip the extension
[~, nameOnly, ~] = fileparts(table_file_name);

% Force the output to be .xlsx
cleanName = regexprep(nameOnly, '\.', '_'); 
cleanName = strtrim(cleanName);
outfile_name = ['PROCESSED_', cleanName, '.xlsx'];

outFile = fullfile(outDir, outfile_name);

%% File Statistics
% Get all IDs of cells contained in the file (should be 1-10)
uniqueIDs = unique(data.TID);
% Make an iterable based on length of uniqueIDs
numTracks = length(uniqueIDs);

% Determine maximum track length for building the MSD table
trackCounts = groupsummary(data, 'TID');
maxLags = max(trackCounts.GroupCount) - 1;

% WARNING: ASSUMES CONSISTENT TIME STEP INTERVALS
timeStep = data.t_min_(2) - data.t_min_(1);

%% Output Setup
%  Set up an empty table for storing speed, velocity, persistence
stats = table( ...
        nan(numTracks,1), nan(numTracks,1), nan(numTracks,1), ...
        nan(numTracks,1), nan(numTracks,1), nan(numTracks,1), ...
        nan(numTracks,1), nan(numTracks,1), nan(numTracks,1), ...
        nan(numTracks,1), nan(numTracks,1), strings(numTracks,1), ...
        'VariableNames', { ...
        'Cell_ID', 'O_Speed', 'O_Velocity', ...
        'O_Persistence', 'O_Displacement', 'N_Speed', ...
        'N_Persistence', 'N_SE', 'N_GoF', ...
        'Fitted_Gamma', 'Fitted_Alpha','PRW_Fit_Flag' ...
        });

% Set up an empty table for storing MSD of each cell
msdVarNames = [{'Tau_min'}, arrayfun(@(x) sprintf('TID_%d', x), uniqueIDs', 'UniformOutput', false)];
msd_table = array2table(nan(maxLags, numTracks + 1), 'VariableNames', msdVarNames);
% Fill in all possible tau values
msd_table.Tau_min = (1:maxLags)' * timeStep;

%% Main Loop
for i = 1:numTracks

    disp(['Processing cell: ', num2str(i)])

    % Downselect to the current cell
    currentTrack = data(data.TID == uniqueIDs(i), :);

    % If length of track is less than min_length, then skip analysis
    if height(currentTrack) <= min_length
        % All values were initialized as NaN, so we just continue
        stats.Cell_ID(i)        = i;
        stats.PRW_Fit_Flag(i)       = "No_Fit";
        continue
    end
    
    % Get the total distance traveled (micron)
    totalDistance = currentTrack.Len_micron_(end);

    % Calculate the elapsed time of cell track (all in minutes)
    timeStart       = currentTrack.t_min_(1);
    timeEnd         = currentTrack.t_min_(end);
    delta_t         = currentTrack.t_min_(3) - currentTrack.t_min_(2) ;
    elapsedTime     = timeEnd - timeStart;

    % Determine start and end points (all in microns)
    startX      = currentTrack.x_micron_(1);
    startY      = currentTrack.y_micron_(1);
    endX        = currentTrack.x_micron_(end);
    endY        = currentTrack.y_micron_(end);

    % Calculate final total displacement (in microns)
    displacement = sqrt((endX - startX)^2 + (endY - startY)^2);
    
    % Calculate MSD and populate table (microns^2)
    cellMSD = calculateMSD(currentTrack.x_micron_, currentTrack.y_micron_);
    msd_table{1:length(cellMSD), i+1}  = cellMSD;

    % Model fitting
    taus = (1:length(cellMSD))' * timeStep;

    % Fit PRW using MSD data
    [fit_persistence,fit_speed,SE,gof] = msd2pse0(cellMSD, delta_t, 2);
    r_squared = gof.rsquare;

    % Fit alpha and Gamma parameters for power law
    [alpha, gamma] = fitAlpha(taus, cellMSD);

    % Append values to statistical table
    stats.Cell_ID(i)          = i;  % Integer

    % Legacy Calculation Methods
    stats.O_Speed(i)          = totalDistance / elapsedTime * 60; % um/hr
    stats.O_Velocity(i)       = displacement / elapsedTime * 60; % um/hr
    % This persistence produces a dimensionless quantity, some ratio
    stats.O_Persistence(i)    = stats.O_Velocity(i) / stats.O_Speed(i); % -
    stats.O_Displacement(i)   = displacement;   % um

    % Model Fit Parameters
    stats.Fitted_Gamma(i)     = gamma;
    stats.Fitted_Alpha(i)     = alpha;
    stats.N_Speed(i)          = fit_speed * 60;   % um / hr
    % Whereas this persistence is a calculated quantity with units. Here,
    % it is converted and report in hours
    stats.N_Persistence(i)    = fit_persistence / 60;   % hr
    stats.N_SE(i)             = SE;
    stats.N_GoF(i)            = r_squared;      % R^2
    if r_squared >= target_r2_val
        stats.PRW_Fit_Flag(i)     = "Pass";
    else
        stats.PRW_Fit_Flag(i)     = "Fail";
    end

end

writetable(stats, outFile, 'Sheet', 'ParamSummary')
writetable(msd_table, outFile, 'Sheet', 'MSD_Data')
disp(stats)

toc

end

%% -- LOCAL FUNCTIONS --

function [alpha, gamma] = fitAlpha(taus, cellMSD)
    % FITALPHA fits the power law: MSD = Gamma * tau^alpha
    % Inputs: taus (vector), cellMSD (vector)
    % Output: alpha (exponent), gamma (anamolous diffusivity parameter)

    validIdx = ~isnan(cellMSD) & cellMSD > 0;
    x = taus(validIdx);
    y = cellMSD(validIdx);

    if length(x) <3
        alpha = NaN;
        gamma = NaN;
        return;
    end

    % Generate intial guess values from taking log of the power law
    p = polyfit(log(x), log(y), 1);
    initial_alpha = p(1);
    initial_gamma = exp(p(2));

    % Fit power law
    ft = fittype('c * x^a', 'independent', 'x', 'dependent', 'y');
    opts = fitoptions(ft);
    opts.StartPoint = [initial_gamma, initial_alpha];
    opts.Lower = [0.001, 1e-10];
    opts.Upper = [2, 1e4];
    opts.Display = 'off';

    try
        fitresult = fit(x(:), y(:), ft, opts);
        alpha = fitresult.a;
        gamma = fitresult.c;
    catch
        alpha = NaN;
        gamma = NaN;
    end
end

function msd = calculateMSD(x, y)
    % Calculates the Mean Squared Displacement for a single track
    % x, y: vectors of coordinates

    N = length(x);
    msd = nan(N-1, 1);

    for tau = 1:(N-1)
        % Calculate displacement for lag tau
        % (Position(t+tau) - Position(t)
        dx = x(1+tau:end) - x(1:end-tau);
        dy = y(1+tau:end) - y(1:end-tau);

        squared_displacements = dx.^2 + dy.^2;
        msd(tau) = mean(squared_displacements, 'omitnan');
    end
end

function [P,S,SE,gof]=msd2pse0(msd,dt,dim)
%******************************************************
% Characterizing the msd by persistent Random walk model.
% output: persistence, speed, and positioning error, goodness of fit.
% input: msd (unit length)
%        dt ( time step (time/frame))
%*******************************************************
% written by : 
%               Pei-Hsun Wu,   PhD
%               Department of Chemical and biomolecular Engineering
%               Johns Hopkins University
% 
% Last update:  Mar, 22, 2012
%               
%*****************************************************
%
% ~~~~~~~~~~THIS FUNCTION HAS BEEN MODIFIED~~~~~~~~~~~~~ 
% List of modifications includes:
%   - Modified upper and lower thresholds
%   - Inclusion of coefficient labels to clarify parameter bound assignments
%   - Truncated data series to be used for model fitting is the first HALF 
%   of the original data series 

    if nargin==1
        dt=1;
    end
    dispfres=0;
    mout=msd;
    toi=1:round(length(mout)/2);
    ti=(1:length(mout(:,1))) ;
    ti0=ti;
    ti0=ti0*dt;
    ti=ti(:);
    ti=ti(toi);

    Nt=length(ti);
    wif=(2*ti.^2+1)/3./ti./(Nt-ti+1) ;    % resolution of MSD at different time lag
    ti=ti.* dt; % convert to real unit

    s = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0 0.01 0],...
                   'Upper',[3 60 100],...
                   'Startpoint',[1 10 1]);
               
    f=fittype([num2str(dim),'*S^2*P*(x-P*(1-exp(-x/P)))+',num2str(2*dim),'*se'],...
        'options',s, ...
        'coefficients', {'S', 'P', 'se'});    

    MSD=mout(toi);
    wt=1./wif.^2./ MSD ; 
    [c2,gof]=fit(ti,MSD,f,'weight',wt) ;            
    yfit=c2(ti);

    % display fitting result 
    if dispfres==1
        figure(999);
        loglog(ti,MSD,'b.',ti,yfit,'r-','linewidth',2); bjff3
        hold on;
        loglog(ti0,mout,'b.',ti0,c2(ti0),'r-','linewidth',2); bjff3
        hold off;
    end
        
    P=c2.P;
    S=c2.S;
    SE=c2.se;
    
end
