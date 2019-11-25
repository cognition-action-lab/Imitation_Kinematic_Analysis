%script to update and merge marks into a compiled data block

homepath = '~/Desktop/';


fprintf('\nSelect new Compiled Data file.\n');
[fnamecompiled,fpathcompiled] = uigetfile('*.*','Select Compiled data file');

fprintf('\nSelect old Marked Data file.\n');
[fnamemarked,fpathmarked] = uigetfile('*.*','Select Marked data file');

tmpcompiled = load(fullfile(fpathcompiled,fnamecompiled));

tmpmarked = load(fullfile(fpathmarked,fnamemarked));

%% post-process the compiled data file to match the marked data

%do some additional pre-processing to the compiled data
for c = 1:size(tmpcompiled.BlockData.Data,1)
    for c2 = 1:size(tmpcompiled.BlockData.Data,2)
        tmppos = [];
        tmprotang = [];
        tmpang = [];
        for d = 1:length(tmpcompiled.BlockData.Data(c,c2).m)
            tmppos = cat(3,tmppos,tmpcompiled.BlockData.Data(c,c2).m(d).pos);
            tmprotang = cat(4,tmprotang,tmpcompiled.BlockData.Data(c,c2).m(d).rotang);
            tmpang = cat(3,tmpang,tmpcompiled.BlockData.Data(c,c2).m(d).angs);
        end
        tmpcompiled.BlockData.pos{c,c2} = tmppos; %*scalefactor;
        tmpcompiled.BlockData.ja{c,c2} = BlockData.Data(c,c2).ja;
        tmpcompiled.BlockData.rotang{c,c2} = tmprotang;
        tmpcompiled.BlockData.ang{c,c2} = tmpang;
        if isfield(tmpcompiled.BlockData.Data(c,c2),'time')
            tmpcompiled.BlockData.time{c,c2} = tmpcompiled.BlockData.Data(c,c2).time;
        elseif isfield(tmpcompiled.BlockData.Data(c,c2),'TrackerTime')
            tmpcompiled.BlockData.time{c,c2} = tmpcompiled.BlockData.Data(c,c2).TrackerTime;
        elseif isfield(tmpcompiled.BlockData.Data(c,c2),'Time')
            tmpcompiled.BlockData.time{c,c2} = tmpcompiled.BlockData.Data(c,c2).Time;
        end
    end
end

%remove the bulky Data field now that we've extracted pos and angle info
tmpcompiled.BlockData = rmfield(tmpcompiled.BlockData,'Data');

%add velocity to dataset
tmpcompiled.BlockData = addvelsubj(tmpcompiled.BlockData);


%% transfer the marks

%copy over the marked indices
tmpcompiled.BlockData.inds = tmpmarked.BlockData.inds;


%% save the output

savfname = fnamemarked;
if ~input('Replace existing marked data file? (Y = 1, N = 0): ')
    savfname = strrep(savfname,'markdata','markdata1'); %create a unique filename to save to, otherwise we will overwrite the existing file
end

BlockData = tmpcompiled.BlockData;

save(fullfile(fpathmarked, savfname),'BlockData');

        
