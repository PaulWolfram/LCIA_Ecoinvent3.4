%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script builds a matrix-based LCA model using the Ecoinvent 3.4 database 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% By Paul Wolfram
% Last updated 22 Aug 2018

clc
clear

%% Step 1: Create meta files
[num txt raw] = xlsread('Inputs/Hybrid_USEEIO_EXIO/EI34_ALabs.xlsx'); % Ecoinvent matrix labels (renamed from the original file "ie_index")
EI_meta.ALabs = [raw(2:14928,6), txt(2:14928,5)];

[num txt raw] = xlsread('Inputs/Hybrid_USEEIO_EXIO/EI34_CLabs.xlsx'); % characterization matrix labels (renamed from the original file "LCIA_index")
EI_meta.CLabs = txt(2:747,4);

[num txt raw] = xlsread('Inputs/Hybrid_USEEIO_EXIO/EI34_SLabs.xlsx'); % normalized emissions matrix labels (renamed from the original file "ee_index")
EI_meta.SLabs = txt(2:2029,9);

clear num raw txt


%% Step 2: Load technology matrix A
EI.A = csvread('Inputs/Hybrid_USEEIO_EXIO/EI34_A.csv'); % load technology matrix from Ecoinvent (renamed from the original file "A_public")
    % The original file "A_public" (provided by Ecoinvent) has been manipulated in Excel beforehand:
    % 1. Text per cell has been split to several columns
    % 2. Everything other than row and column indices and the coefficients
    % has been deleted
    % 3. Indices were changed to start at 1 instead of at 0
EI.A(:,1:3) = EI.A;                   % only keep first 3 columns
EI.A = spconvert(EI.A);               % convert sparse list to sparse matrix 

EI.A = full(EI.A);
EI.A(isnan(EI.A))=0;
EI.A(isinf(EI.A))=0;

% According to readme file, off-diagonal elements have opposite sign, so
EI.A = -EI.A;                         % ... so I revert all signs ...           
EI.A(1:(length(EI.A)+1):end) = 0; % and set diagonal elements to zero


%% Step 3: Load characterization matrix C
EI.C = csvread('Inputs/Hybrid_USEEIO_EXIO/EI34_C.csv'); % load technology matrix from Ecoinvent (renamed from the original file "C")
EI.C(:,1) = EI.C(:,1) + ones(size(EI.C,1),1); % change first process number to 1 (EI processes start at number 0)
EI.C(:,2) = EI.C(:,2) + ones(size(EI.C,1),1); 
EI.C = spconvert(EI.C); % convert to sparse matrix 
EI.C = full(EI.C);


%% Step 4: Load emissions factors S
EI.S = csvread('Inputs/Hybrid_USEEIO_EXIO/EI34_S.csv'); % load technology matrix from Ecoinvent (renamed from the original file "B_public")
EI.S(:,1) = EI.S(:,1) + ones(size(EI.S,1),1); % change first process number to 1 (EI processes start at number 0)
EI.S(:,2) = EI.S(:,2) + ones(size(EI.S,1),1); 
EI.S = EI.S(:,1:3); % keep only first 3 cols
EI.S = spconvert(EI.S); % convert to sparse matrix 
EI.S = full(EI.S);


%% Step 5: add some meta data
EI_meta.Nprocesses =  size(EI.A,1); % number of processes
EI_meta.Nemissions =  size(EI.S,1); % number of emissions
EI_meta.Nimpact_cats =  size(EI.C,1); % number of impact categories


%% Step 6: Convert all GHG emissions to CO2e
a = strmatch('IPCC 2013;climate change;GWP 100a;kg CO2-Eq', EI_meta.CLabs(:,1)); % identify GWP 100
EI.S_GHG = EI.C(a,:) * EI.S; 
EI.S_GHG = full(EI.S_GHG)
EI.S_GHG = [EI.S_GHG, zeros(1, size(EI.A,1)-size(EI.S_GHG,2))] % fill the last empty elements
 

%% Step 7: Create a final demand vector y
EI.y = zeros(size(EI.A,1),1);


%% Step 8: Calculate Leontief inverse
% This step takes a few minutes
L = inv(eye(size(EI.A)) - EI.A); % calculate Leontief inverse


%% Step 8: Example 1: Coal electricity in Austria
a = strmatch('electricity production  hard coal; AT; kWh', EI_meta.ALabs(:,2)); % identify the process, compare with ALabs
y_coal = EI.y;
y_coal(a) = 1; % define a demand of 1 (kWh)
CF_coal =  EI.S_GHG * L * y_coal; % calculate carbon footprint
process_name = EI_meta.ALabs(a,2);
fprintf('"%s" causes direct industry GHG emissions of %.3f kg CO2e/kWh.\n', process_name{1}, EI.S_GHG(a)); % Print answer
fprintf('"%s" causes life-cycle GHG emissions of %.3f kg CO2e/kWh.\n', process_name{1}, CF_coal); % Print answer

% Contribution analysis:
CF_coal_contr =  diag(EI.S_GHG) * L * y_coal; % repeat CF calc w/ diagonalized S vector
[value index] = sort(CF_coal_contr,'descend'); % sort biggest contributors in descending order
top25_values = value(1:25) % get top 25 
top25_indices = index(1:25) % get their indices
top25_sectors = EI_meta.ALabs(top25_indices,2) % get process names


%% Step 9: Example 2: BEV, global average
clc
clear a CF_coal y_coal CF_coal_contr
a = strmatch('transport  passenger car  electric; GLO; km', EI_meta.ALabs(:,2)); % identify the process
y_BEV = EI.y;
y_BEV(a) = 1;
CF_BEV =  EI.S_GHG * L * y_BEV; % calculate carbon footprint
process_name = EI_meta.ALabs(a,2);
fprintf('"%s" causes direct industry GHG emissions of %.3f kg CO2e/km.\n', process_name{1}, EI.S_GHG(a)); % Print answer
fprintf('"%s" causes life-cycle GHG emissions of %.3f kg CO2e/km.\n', process_name{1}, CF_BEV); % Print answer

% Contribution analysis:
CF_BEV_contr =  diag(EI.S_GHG) * L * y_BEV; % repeat CF calc w/ diagonalized S vector
[value index] = sort(CF_BEV_contr,'descend'); % sort biggest contributors in descending order
top25_values = value(1:25) % get top 25 
top25_indices = index(1:25) % get their indices
top25_sectors = EI_meta.ALabs(top25_indices,2) % get process names


%% Step 10: Example 3: ICEV, European average
clc
clear a CF_BEV y_BEV CF_BEV_contr
a = strmatch('transport  passenger car with internal combustion engine', EI_meta.ALabs(:,2)); % identify the process
y_ICEV = EI.y;
y_ICEV(a(1)) = 1;
CF_ICEV =  EI.S_GHG * L * y_ICEV; % calculate carbon footprint
process_name = EI_meta.ALabs(a,2);
fprintf('"%s" causes direct industry GHG emissions of %.3f kg CO2e/km.\n', process_name{1}, EI.S_GHG(a(1))); % Print answer
fprintf('"%s" causes life-cycle GHG emissions of %.3f kg CO2e/km.\n', process_name{1}, CF_ICEV); % Print answer

% Contribution analysis:
CF_ICEV_contr =  diag(EI.S_GHG) * L * y_ICEV; % repeat CF calc w/ diagonalized S vector
[value index] = sort(CF_ICEV_contr,'descend'); % sort biggest contributors in descending order
top25_values = value(1:25) % get top 25 
top25_indices = index(1:25) % get their indices
top25_sectors = EI_meta.ALabs(top25_indices,2) % get process names


%% Step 11: Example 4: Market for ICEV, small, EURO5, European average
clc
clear a CF_ICEV y_ICEV CF_ICEV_contr
a = strmatch('market for transport  passenger car  small size  petrol  EURO 5', EI_meta.ALabs(:,2)); % identify the process
y_ICEVmkt = EI.y;
y_ICEVmkt(a(1)) = 1;
CF_ICEVmkt =  EI.S_GHG * L * y_ICEVmkt; % calculate carbon footprint
process_name = EI_meta.ALabs(a,2);
fprintf('"%s" causes direct industry GHG emissions of %.3f kg CO2e/km.\n', process_name{1}, EI.S_GHG(a(1))); % Print answer
fprintf('"%s" causes life-cycle GHG emissions of %.3f kg CO2e/km.\n', process_name{1}, CF_ICEVmkt); % Print answer

% Contribution analysis:
CF_ICEVmkt_contr =  diag(EI.S_GHG) * L * y_ICEVmkt; % repeat CF calc w/ diagonalized S vector

[value index] = sort(CF_ICEVmkt_contr,'descend'); % sort biggest contributors in descending order
top25_values = value(1:25) % get top 25 
top25_indices = index(1:25) % get their indices
top25_sectors = EI_meta.ALabs(top25_indices,2)


%% Step 12: Save
clear a ans process_name CF_ICEVmkt CF_ICEVmkt_contr y_ICEVmkt top25_indices top25_sectors top25_values value index
save('EI34.mat','EI','EI_meta') % Note that L is too large to save
disp('Done.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
