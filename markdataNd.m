% Place markers in data, especially saccades.
% Plot data and fast phase markers,
% and interactively edit the fast phase markers.  Three traces are
% displayed, with marks from the changeable (middle) trace also appearing
% on the other two traces for alignment purposes.
%
% ind=markdata3d(x,y,z,sssize,peakflag,ind0,gridflag,Fs,'arg',param)
% x=data
% sssize=screen size - number of points to plot per screen
% peakflag=1: program forces indices to peaks
% ind0 = y; varargin: indx = x, indz = z,
% extra parameters: 'Name', 'n1', 'idx', 'idz'
% note that changes can only be made to the y data set markings (ind0)
% When the data and markers are displayed:
% Left mouse button = add a fast phase marker
% Right mouse button = remove a fast phase marker
% 'n' = advance to next data section
% 'x' = end
% 'r' = redraw screen (needed after several removals)
% 'g' = go to (prompts for sample number to go to)
% 's' = change ssize
% 'd' = find duration and slope of delineated segment nearest cursor
% 't' = find time interval
% '1' = invert upper pannel
% '2' = invert middle pannel
% '3' = invert lower pannel
% '4' = increase y limits of upper pannel
% 'u' = undo remove mark
%

% MS  May 1993
% Jan 1995 - changed button code for deletion from 2 to 3
% Jan 1997 - changes for MATLAB 5 (clf and colors)
% Jan 2007 - changed to single data set

%we will assume that for matrix X, dimension 1 is time, dimension 2 is
%across plots, and dimension 3 (if > 1) is within plots

function ind=markdataNd(X,sssize,peakflag,ind0,gridflag,varargin)

n1 = 1;
resetylim = 0;

for a = 1:2:length(varargin)
    switch(lower(varargin{a}))
        case 'name'
            set(gcf,'Name',varargin{a+1})
        case 'n1'
            n1 = round(varargin{a+1}/1000)*1000+1;
    end
end

ssize=sssize;
N=size(X,1);
Nplots = size(X,2);

%n1=1;
n2=min(N,ssize)+n1-1;

if size(ind0,1) > size(ind0,2)
    ind0 = ind0';
end

ii=ind0;

%figure

% plot and edit until 'x' is pressed
but='n';
while(but~='x')
    % plot position data with markers

    if(but=='n')
        clf;hold off;
        for iplots = 1:Nplots
            h = subplot(Nplots,1,iplots);
            plot((n1:n2),squeeze(X(n1:n2,iplots,:)),'-')
            ii1=ii(1:2:length(ii));
            ii2=ii(2:2:length(ii));
            iii1=ii1(find(ii1>=n1 & ii1<=n2));
            iii2=ii2(find(ii2>=n1 & ii2<=n2));
            plot((n1:n2),squeeze(X(n1:n2,iplots,:)),'-')
            hold on
            for ilines = 1:size(X,3)
                if(~isempty(iii1)) plot(iii1,X(iii1,iplots,ilines),'xr'); end
                if(~isempty(iii2)) plot(iii2,X(iii2,iplots,ilines),'or'); end
            end
            if(gridflag) grid on; end
        end            
        hold off;
        if resetylim == 1
            ylims = get(h,'Ylim');
            ylims = ylims*1.1;
            set(h,'Ylim',ylims);
        end
        
    end

    % edit the indices
    [xx,yy,but]=ginput(1);
    % left button - add a point
    if(but==1)
        inew=fix(xx);
        if(inew>5 & inew<N-5)
            aa=y(inew)-y(inew-2);
            yseg=y(inew-3:inew+3);
            if(peakflag==1)
                if(aa>=0) inew=inew-4+find(yseg==max(yseg)); inew=inew(1); end
                if(aa<0)  inew=inew-4+find(yseg==min(yseg)); inew=inew(1); end
            end
        end
        ii=[ii inew];
        ii=sort(ii);
        for iplots = 1:Nplots
            subplot(Nplots,1,iplots);
            for ilines = 1:size(X,3)
                hold on; 
                plot(inew,X(inew,iplots,ilines),'or'); 
                hold off;
            end
        end  
        
    end

    % right button - delete a point
    if(but==3)
        dii=abs(ii-fix(xx));       % distances from cursor to each index
        iomit=find(dii==min(dii)); % min dist is the one to remove
        omit=ii(iomit);
        ii=[ii(1:iomit-1) ii(iomit+1:length(ii))];
        for iplots = 1:Nplots
            subplot(Nplots,1,iplots);
            for ilines = 1:size(X,3)
                hold on; 
                plot(omit,X(omit,iplots,ilines),'xw',omit,X(omit,iplots,ilines),'ow'); 
                hold off;
            end
        end  
        
    end

    % key 'n' - advance to next data section
    if(but=='n')
        n1=n2+1;
        if(n1>=N) n1=N-ssize; end
        n2=min(n1+ssize-1,N);
    end

    % key 'b' - back to last data section
    if(but=='b') n1=max(1,n1-ssize); n2=min(n1+ssize-1,N); but='n'; end

    % key 'r' - redraw screen
    % act as if 'n' pressed, but don't change indices
    if(but=='r') but='n'; end
    
    %key 'u' - undo remove mark
    if(but=='u')
        if isempty(find(ii==omit))
            ii = sort([ii omit]);
        end
        but = 'n';
    end        
    
    % key 's' - change ssize
    % then act as if 'n' pressed, but don't change indices
    if(but=='s')
        ssize=input('New value for ssize: ');
        n2=min(n1+ssize-1,N);
        but='n';
    end

    % key 'g' - go to specified section
    if(but=='g')
        n1=input('Go to sample number: ');
        if(n1>=N) n1=N-ssize; end
        n2=min(n1+ssize-1,N);
        but='n';
    end

    
end

hold off
subplot
ind=ii;
return
