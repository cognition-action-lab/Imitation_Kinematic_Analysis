%batch process script for ABblockcompile


paths = {
            %**********************************************************************************
            %** leave this empty unless processing files located in different parent folders **
            %**   in that case, hard code the paths below.                                   **
            %**********************************************************************************
            
            %'Data/Named Meaningful/Controls/41/'
            %'Data/Named Meaningful/Controls/116/'
            %'Data/Named Meaningful/Controls/344/'
            %'Data/Named Meaningful/Controls/345/'
           };

       
if isempty(paths)
    fprintf('\nSelect all data folders to be analyzed.\n');
    paths = uigetmultidir();  %call GUI to allow multiselect folders
end

if isempty(paths)
    fprintf('\n\nERROR: No data files selected for analysis.\n\n');
    return;
end

[savepath] = uigetdir('','Choose folder to save individual subject data into (cancel to save into the individual subject folder)');
if isempty(savepath)
    savepath = 0;
end

% fprintf('\nSelect data file to append to (or cancel to create a new data file).\n');
% [datafname,datafpath] = uigetfile('*.mat','Select data file to append to.');


for a = 1:length(paths)
    ABblockcompile(paths{a},savepath);
end

%ABsubjectcompile(paths,fullfile(datafpath,datafname));
