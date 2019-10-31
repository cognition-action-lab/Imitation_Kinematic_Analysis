
clear all;
close all;


fprintf('\nSelect data files to process.\n');
[fnames,fpath] = uigetfile('*.mat','Select data file to process.','Multiselect','on');


for a = 1:length(fnames)
    clear BlockData;
    
    load(fullfile(fpath,fnames{a}));
    
    
    %do some pre-processing since we've skipped the "ABsubjectcompile" code
    
    for b = 1:length(BlockData)
        if isempty(BlockData{b})
            continue;
        end
        
        if b == 10
            scalefactor = 1; %note, this is to convert meters to inches (1/0.0254)
            scaletype = 'scalemeters';
        else
            scalefactor = (0.0254/1);  %note, this is to convert inches to meters (0.0254/1)
            scaletype = 'scaleinches';
        end
        
        for c = 1:size(BlockData{b}.Data,1)
            for c2 = 1:size(BlockData{b}.Data,2)
                tmppos = [];
                tmprotang = [];
                tmpang = [];
                for d = 1:length(BlockData{b}.Data(c,c2).m)
                    tmppos = cat(3,tmppos,BlockData{b}.Data(c,c2).m(d).pos);
                    tmprotang = cat(4,tmprotang,BlockData{b}.Data(c,c2).m(d).rotang);
                    tmpang = cat(3,tmpang,BlockData{b}.Data(c,c2).m(d).angs);
                end
                BlockData{b}.pos{c,c2} = tmppos*scalefactor;
                BlockData{b}.ja{c,c2} = BlockData{b}.Data(c,c2).ja;
                BlockData{b}.rotang{c,c2} = tmprotang;
                BlockData{b}.ang{c,c2} = tmpang;
                if isfield(BlockData{b}.Data(c,c2),'time')
                    BlockData{b}.time{c,c2} = BlockData{b}.Data(c,c2).time;
                elseif isfield(BlockData{b}.Data(c,c2),'TrackerTime')
                    BlockData{b}.time{c,c2} = BlockData{b}.Data(c,c2).TrackerTime;
                elseif isfield(BlockData{b}.Data(c,c2),'Time')
                    BlockData{b}.time{c,c2} = BlockData{b}.Data(c,c2).Time;
                end
            end
        end
        
        %remove the bulky Data field now that we've extracted pos and angle info
        BlockData{b} = rmfield(BlockData{b},'Data');
        
        %add velocity to dataset
        BlockData{b} = addvelsubj(BlockData{b});
        
    end
    
    tmp = load(fullfile(strrep(fpath,'analyzed data2','analyzed data1'),strrep(fnames{a},'_data','_markdata')));
    
    for b = 1:length(tmp.BlockData)
        if isempty(tmp.BlockData{b})
            continue;
        end
        
        BlockData{b}.inds = tmp.BlockData{b}.inds(1:41,:);
    end
    
    save(fullfile(fpath,strrep(fnames{a},'_data','_markdata')),'BlockData');
    
    clear tmp*
    
end



