function multiplot2D(ydata,ydataname,varargin)

% multiplot2D.m  (c) T.Kam 2006
%
% This function plots multiple data series in subplots typically suitable
% for journal layout sizes
%
% data is T x n matrix, where T is #points, n is #variables
% dataname is n x 1 string vector
% varargin is "xdata" for x-axis labels

nwin = size(ydata,2); % #subplots, nwin = n
yscale = 0.2; % create buffer with size of +/- 0 <= yscale <= 1

figure
    if nwin <= 6
        wc = 2;                           % window cols
    elseif nwin > 6 & nwin <= 9
        wc = 3;
    else
        wc = 4;
    end
    
        wr = ceil(nwin/wc); % window rows
    
    for i = 1:nwin     
        subplot(wr,wc,i)
            
            if nargin <= 2
                plot(ydata(:,i));
            elseif nargin > 2
                plot(varargin{:},ydata(:,i));
            end
            ylabel(ydataname(i),'FontSize',16)
        grid on
       
        axis([min(varargin{:}), max(varargin{:}), ...
              (min(ydata(:,i)) - (max(ydata(:,i))-min(ydata(:,i)))*yscale), ...
              (max(ydata(:,i))+(max(ydata(:,i))-min(ydata(:,i)))*yscale)])
    end   