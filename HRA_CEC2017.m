clc;        % Clear Command Window
clear;      % Clear workspace variables
close all;  % Close all figures

% Check if 'HRA_CEC2017.txt' already exists
outfile ='HRA_CEC2017.txt'; 
if exist(outfile, 'file') == 2
    delete(outfile); % Remove the file if it exists
end
% Start capturing the command window output
diary(outfile);
diary on;
% Define the Excel output file name
excelFileName = fullfile('HRA_CEC2017_results.xlsx');


% Define the folder containing the Excel files
folderName = 'CEC2017_Results';
fileNames = {'CEC2017_Dim_10.xlsx', 'CEC2017_Dim_30.xlsx', 'CEC2017_Dim_50.xlsx', 'CEC2017_Dim_100.xlsx'};
numDimensions = length(fileNames);

% Define the names of the alternatives
SelAlgorithms = {...
    'jSO', 'MM_OED', 'IDEbestNsize', 'RB_IPOP_CMA_ES', 'LSHADE_SPACMA', 'DES', ...
    'DYYPO', 'TLBO_FL', 'PPSO', 'MOS_SOCO2011', 'MOS_CEC2013', 'LSHADE_cnEpSin', 'EBOwithCMAR'};

NoAlgorithms = length(SelAlgorithms);
% Initialize matrices to store ranks
rankMatrix_D = zeros(NoAlgorithms, numDimensions); % 13 algorithms


% Initialize matrices to store ranks
rankMatrix_D = zeros(13, numDimensions); % 13 algorithms
decisionMatrices = cell(1, numDimensions); % Cell array to store decision matrices for each dimension

for d = 1:numDimensions
    % Get the file name without extension (dimension name)
    [~, fileName, ~] = fileparts(fileNames{d});

    % Construct the full file path
    filePath = fullfile(folderName, fileNames{d});

    % Get the sheet names from the Excel file
    [~, sheetNames] = xlsfinfo(filePath);

    % Initialize matrix for performance measures ranks
    ranks_P = zeros(13, 5); % 13 algorithms, 5 performance measures

    for p = 1:5
        % Read the data from the p-th sheet
        sheetData = readmatrix(filePath, 'Sheet', p);

        % Transpose the sheetData to have algorithms as rows and functions as columns
        sheetData = sheetData';

        % Initialize matrix to store ranks
        rankMatrix = zeros(size(sheetData));

         % Define criteria signs (+1 for benefit, -1 for cost)
        criteriaSign = -ones(1, size(sheetData, 2)); % All criteria are costs

        % Compute ranks for each column
        for col = 1:size(sheetData, 2)
            rankMatrix(:, col) = tiedrank(sheetData(:, col));
        end
        % set the input ordinal matrix 
        N=rankMatrix;
        %disp(N);
        
        [m,n]=size(N);
        criteriaSign = (-1)*ones(1,n); % all criteria represent cost
        % Choose the normalization method for R-TOPSIS
        method = 'Max-Min';%'Max'; % Choose between 'Max' or 'Max-Min'
        % fixed domain [0,m+1] : Max-Min & Max are equivalent 
        domain = zeros(2,n);
        domain(2,:)= (m+1)*ones(1,n); % upper bound is  m=#Algorithms
        weights_F = (1/n)*ones(1,n); 

        % Call the RTOPSIS function
        [sortedScore, ~, sortedIndices, Rank, ~] = ...
            RTOPSIS(N, weights_F, criteriaSign, domain, method, SelAlgorithms);

        fprintf('=============================================================\n')
        disp(['R-TOPSIS results for dimension: ' fileName ', Sheet: ' sheetNames{p} ':']);
        %disp(RTOPSIS_tbl)
        
        % Display the ranking for the original list of algorithms
        % Initialize the original rankings and scores arrays
        OriginalRanking = zeros(size(SelAlgorithms));
        OriginalScores = zeros(size(SelAlgorithms));

        % Use the sorted indices to map the ranks and scores back to the original order
        OriginalRanking(sortedIndices) = Rank;
        OriginalScores(sortedIndices) = sortedScore;
   
        % Display Original Rankings and Scores
        T = table(SelAlgorithms(:), OriginalScores(:), OriginalRanking(:), 'VariableNames', {'Algorithm', 'Closeness', 'Rank'});
        disp(T);
        
        % Store the ranks for this performance measure
        ranks_P(:, p) = OriginalRanking;

        % Store the ranks in the decision matrix
        decisionMatrix(:, p) = OriginalRanking;

    end

    % Store the decision matrix for this dimension
    decisionMatrices{d} = decisionMatrix;
    fprintf('=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*\n')
    % Display the decision matrix for this dimension
    disp(['Decision matrix for dimension ' fileName ':']);
    disp(decisionMatrix);
    disp(' ');

    % Apply R-TOPSIS to the ordinal decision matrix
    
    [m,n]=size(decisionMatrix);
    criteriaSign = (-1)*ones(1,n); % all criteria represent cost
    % Choose the normalization method for R-TOPSIS
    method = 'Max-Min';%'Max'; % Choose between 'Max' or 'Max-Min'
    % fixed domain [0,m+1] : Max-Min & Max are equivalent
    domain = zeros(2,n);
    domain(2,:)= (m+1)*ones(1,n); % upper bound is  m=#Algorithms
    weights_P = (1/n)*ones(1,n);
    
    % Call the RTOPSIS function
    [sortedScore, ~, sortedIndices, Rank, ~] = ...
        RTOPSIS(decisionMatrix, weights_P, criteriaSign, domain, method, SelAlgorithms);

    fprintf('=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*\n')
    disp(['R-TOPSIS results for dimension: ' fileName ]);
    disp('weights for each performance measure: [best, worst, median, mean, std]')
    disp(weights_P)

    % Display the ranking for the original list of algorithms
    % Initialize the original rankings and scores arrays
    OriginalRanking = zeros(size(SelAlgorithms));
    OriginalScores = zeros(size(SelAlgorithms));

    % Use the sorted indices to map the ranks and scores back to the original order
    OriginalRanking(sortedIndices) = Rank;
    OriginalScores(sortedIndices) = sortedScore;

    % Display Original Rankings and Scores
    T = table(SelAlgorithms(:), OriginalScores(:), OriginalRanking(:), 'VariableNames', {'Algorithm', 'Closeness', 'Rank'});
    disp(T);
    
    % Store the ranks
    rankMatrix_D(:, d) = OriginalRanking;
end

% Final decision matrix (13x4) with columns the ranks of the previous TOPSIS
finalDecisionMatrix = rankMatrix_D;

% Display the final decision matrix
fprintf('\nFinal decision matrix D:\n');
disp(finalDecisionMatrix);
disp(' ');

[m,n]=size(finalDecisionMatrix);
criteriaSign = (-1)*ones(1,n); % all criteria represent cost
% Choose the normalization method for R-TOPSIS
method = 'Max-Min';%'Max'; % Choose between 'Max' or 'Max-Min'
% fixed domain [0,m+1] : Max-Min & Max are equivalent
domain = zeros(2,n);
domain(2,:)= (m+1)*ones(1,n); % upper bound,   m=#Algorithms
weights_D = (1/n)*ones(1,n);

% Call the RTOPSIS function
[sortedScore, sortedAlgorithms, sortedIndices, Rank, RTOPSIS_tbl] = ...
    RTOPSIS(finalDecisionMatrix, weights_D, criteriaSign, domain, method, SelAlgorithms);


fprintf('=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*\n')
disp('R-TOPSIS results  - overall ranking : ');
disp('weights for each dimension d=[10, 30, 50, 100]')
disp(weights_D)


% Display the ranking for the original list of algorithms
% Initialize the original rankings and scores arrays
OriginalRanking = zeros(size(SelAlgorithms));
OriginalScores = zeros(size(SelAlgorithms));

% Use the sorted indices to map the ranks and scores back to the original order
OriginalRanking(sortedIndices) = Rank;
OriginalScores(sortedIndices) = sortedScore;

% Display Original Rankings and Scores
T = table(SelAlgorithms(:), OriginalScores(:), OriginalRanking(:), 'VariableNames', {'Algorithm', 'Closeness', 'Rank'});
disp(T);


%%-------------------------------------------------------------------------
% Stop capturing the command window output
diary off;

% Read the captured text
fileID = fopen(outfile, 'r');
outputText = fread(fileID, '*char')';
fclose(fileID);

% Format the text as HTML
htmlText = ['<html><head><title>Command Window Output</title></head><body>', ...
            '<h1>Output from MATLAB Command Window</h1>', ...
            '<pre>', outputText, '</pre>', ...
            '</body></html>'];
% Get the name of the current script
scriptName = mfilename;

% Construct HTML filename using the script name
htmlFileName = sprintf('%s.html', scriptName);

% Save the HTML text to an HTML file named after the current script
htmlFileID = fopen(htmlFileName, 'w');
fprintf(htmlFileID, '%s', htmlText);
fclose(htmlFileID);

% Open the HTML file automatically in the default web browser
web(htmlFileName, '-browser');

