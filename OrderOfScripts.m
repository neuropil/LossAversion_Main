

1. FOOOF_table_create
% Takes the processed ephys and runs FOOOF on it and creates a table and
% saves it 

2. FOOOF_table_power_v1
% Takes FOOOF table from FOOOF_table_create and creates a table for where
% each power and frequency are. Then concatenates all the tables together
- FOOOF_table_power_Blocks_v1
    % Takes the FOOOF table from FOOOF_table_create and creates a table
    % with power and frequency and keeps the block number. Created an inner
    % function called processBlock.m that does inner processing
    % Saves files as partID_hemi_BA_blocks.mat 
    % saved in Y:\LossAversion\LH_Data\FOOOF_data\CLASE001\Tables
   

3. CombineTables
% takes all the tables from all participants and puts them into one table 
    % For FOOOF_table_power_Blocks_v1, after all participants have been combined 
    % it was saved as FOOOF_allBlockTables.mat in
    % Y:\LossAversion\LH_Data\FOOOF_data. It gets loaded as
    % allPartTablesBlocks. I added LaVal and saved the new table as PeriodicByBlocks in 
    % Y:\LossAversion\LH_Data\FOOOF_data\PeriodicComponents. 
        % PeriodicByBlocks has the periodic components by blocks and fooof
        % was ran on each epoch 




% Figure generation 
FOOOF_scatter_v1
- Creates scatter plot of FOOOF power data on y axis and frequency on x axis 
- Creates box plots for aperiodic components 

FOOOF_scatter_v2
- creates boxplots of periodic data 

Frequency_fig_epochs
- Creates scatter of the output of all of the participant data from combine tables 

Frequency_fig_epochs_blocks
- Creates scatter of the output of all of the participant data with the block numbers 




% Periodic Figure Generation 
PeriodicBoxcharts_v1 
% creates boxcharts for each frequency band and their loss aversion
% grouping


%% APERIODIC 

1. WholeExperimentEphys_v1

2. FOOOF_allEphys

3. aperiodicTAB_v1 

4. CombineTables 

5. Aperiodic_Boxchart

%% PSD Fig 

1. FOOOF_PSD_Create

2. combinePSDs

3. CombineTables

4. FigStart_v1

%% FOOOF on each block 

1. FOOOF_byBlocks

2. FOOOF_Table_Blocks

3. CombineFOOOFTables

%% Baseline 

1. baseline_between_blocks_v1

2. baseline_runFOOOF

3. aperiodicTAB_baseline_v1

4. CombineTables 


%% Zscore aperiodic 
%%Trial by trial 
1. FOOOF_trial_by_trial 
% runs fooof on each trial per brain area 

2. FOOOF_inbetween_blocks
% calculates the mean and standard deviation for each break and then summs
% them 

3. Zscore_FOOOF
% Zscores each trial 

4. CombineTables 


%% Granger Causality 
Scripts are in the GrangerCausality folder 

Kevin ran the initial Granger Causality. I think he ran test_GC to create everything. 

1. test_GC
- Located in GrangerCausality/KevinCode
- I think this is what Kevin used to generate the GC values. 
- Uses EpochEphys.mat file found in ProcessedEphys folder 
- Also need the inner function Autoregressive_Process_V1.mat (located in GrangerCausality main folder)

2. OutDegGC
- Calculates the OutDegree from Centrality for each GC
- Saves outdegree and outdegree table to the same file that was loaded 

3. CombineTables
- Saves all participant tables into one large one 


Make GC Heatmap Figures 
- To make heatmap GC figures, use the function Full_Plot in KevinCode folder 

