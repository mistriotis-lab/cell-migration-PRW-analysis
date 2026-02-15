function [] = multi_file_msd_process()
    %% multi_file_msd_process
    % Takes a folder of cell trajectory data files for post processing.
    % Computes speed and persistence time from the mean squared
    % displacement (MSD) data fitted in a persistent random walk (PRW) 
    % model. Additionally, fits the exponent and anamolous diffusivity 
    % coefficient from the power law definition of MSD. 
    %
    % Lines 332 through 334 may be modified to change boundaries and
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
    
    % Get path to folder of data to be analyzed
    folderPath = uigetdir(pwd, 'Select a folder to be analyzed:'); 
    
    % Lists all files in folderPath
    filesStruct = dir(folderPath);
    % Removes the current_dir and parent_dir objects from filesStruct
    filesStruct = filesStruct(~ismember({filesStruct.name}, {'.', '..'}));
    % Removes any subdirectories
    filesStruct = filesStruct(~[filesStruct.isdir]);
    
    % Loop over all remaining files
    for f = 1:length(filesStruct)
    
        % Holds just the file's name 
        baseFileName = filesStruct(f).name;
    
        fprintf('Now processing %s\n', baseFileName);
    
        % Holds the full path to the file
        full_table_path = fullfile(filesStruct(f).folder, baseFileName);    
    
        % If an output file from previous analysis exists in directory, skip
        if startsWith(baseFileName, 'PROCESSED')
            continue
        end
    
        % Load data from file
        data = readtable(full_table_path);
    
        % Check if time column label matches anticipated format
        % If not, the label is changed within the data structure
        if ismember('t_mins_', data.Properties.VariableNames)
            data.Properties.VariableNames{'t_mins_'} = 't_min_';
        end
    
        %% Pathing
        % path to an output directory within the data directory
        outDir = fullfile(folderPath, 'temp_out');
        % Print the output directory path
        disp(outDir)
        % Check if outDir exists, create folder if it does not
        if ~exist(outDir, 'dir')
            mkdir(outDir)
        end
    
        % Use fileparts to strip the extension of the data file
        [~, nameOnly, ~] = fileparts(baseFileName);
        
        % Force the output calculation file to be .xlsx
        cleanName = regexprep(nameOnly, '\.', '_'); 
        cleanName = strtrim(cleanName);
        outfile_name = ['PROCESSED_', cleanName, '.xlsx'];
    
        outFile = fullfile(outDir, outfile_name);
    
        %% File Statistics
        % Get all IDs of cells contained in the file
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
        for c = 1:numTracks
        
            disp(['Processing cell: ', num2str(c)])
        
            % Downselect to the current cell
            currentTrack = data(data.TID == uniqueIDs(c), :);
        
            % If length of track is less than min_length, then skip analysis
            if height(currentTrack) <= min_length
                % All values were initialized as NaN, so we just continue
                stats.Cell_ID(c)        = c;
                stats.PRW_Fit_Flag(c)       = "No_Fit";
                continue
            end
            
            % Get the total distance traveled (micron)
            totalDistance = currentTrack.Len_micron_(end);
        
            % Calculate the elapsed time of cell track (all in minutes)
            timeStart       = currentTrack.t_min_(1);
            timeEnd         = currentTrack.t_min_(end);
            delta_t         = currentTrack.t_min_(3) - currentTrack.t_min_(2) ;
            elapsedTime     = timeEnd - timeStart;
        
            % Assign start and end points (all in microns)
            startX      = currentTrack.x_micron_(1);
            startY      = currentTrack.y_micron_(1);
            endX        = currentTrack.x_micron_(end);
            endY        = currentTrack.y_micron_(end);
        
            % Calculate final total displacement (in microns)
            displacement = sqrt((endX - startX)^2 + (endY - startY)^2);
            
            % Calculate MSD and populate table (microns^2)
            cellMSD = calculateMSD(currentTrack.x_micron_, currentTrack.y_micron_);
            msd_table{1:length(cellMSD), c+1}  = cellMSD;
        
            % Model fitting
            taus = (1:length(cellMSD))' * timeStep;
       
            % Fit PRW using MSD data
            [fit_persistence,fit_speed,SE,gof] = msd2pse0(cellMSD, delta_t, 2);
            r_squared = gof.rsquare;
    
            % Fit alpha and Gamma parameters for power law
            [alpha, gamma] = fitAlpha(taus, cellMSD);
        
            % Append values to statistical table
            stats.Cell_ID(c)          = c;  % Integer
        
            % Legacy Calculation Methods
            stats.O_Speed(c)          = totalDistance / elapsedTime * 60; % um/hr
            stats.O_Velocity(c)       = displacement / elapsedTime * 60; % um/hr
            % This persistence produces a dimensionless quantity, some ratio
            stats.O_Persistence(c)    = stats.O_Velocity(c) / stats.O_Speed(c); % -
            stats.O_Displacement(c)   = displacement;   % um
        
            % Model Fit Parameters
            stats.Fitted_Gamma(c)     = gamma;
            stats.Fitted_Alpha(c)     = alpha;
            stats.N_Speed(c)          = fit_speed * 60;   % um / hr
            % Whereas this persistence is a calculated quantity with units. Here,
            % it is converted and report in hours
            stats.N_Persistence(c)    = fit_persistence / 60;   % hr
            stats.N_SE(c)             = SE;
            stats.N_GoF(c)            = r_squared;      % R^2
            if r_squared >= target_r2_val
                stats.PRW_Fit_Flag(c)     = "Pass";
            else
                stats.PRW_Fit_Flag(c)     = "Fail";
            end
        
        end
        
        writetable(stats, outFile,  'Sheet', 'ParamSummary')
        writetable(msd_table, outFile, 'Sheet', 'MSD_Data')
        disp(stats)
    
    end
    
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
% This function is adapted from PH. Wu, A. Giri, and D. Wirtz, 
% "Statistical analysis of cell migration in 3D using the anisotropic 
% persistent random walk model," *Nat Protoc*, vol. 10, no. 3, pp. 517–527,
% Mar. 2015, doi: [10.1038/nprot.2015.030](10.1038/nprot.2015.030)" and is
% used under the MIT License.
%
% MIT License
% 
% Copyright (c) 2026 DeepBioVision
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

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
%          THIS FUNCTION HAS BEEN MODIFIED 
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
                   'Upper',[3.5 20 100],...
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
