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

for a = 1:length(paths)
    ABblockcompile(paths{a});
end

ABsubjectcompile(paths);
