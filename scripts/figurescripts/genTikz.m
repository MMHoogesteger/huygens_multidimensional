function  genTikz(saveName)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
matlab2tikz([saveName '.tex'],'parseStrings',false,...
    'height','\figureheight',...
    'width','\figurewidth',...
    'showInfo', false,...
    'extraAxisOptions', {'scaled ticks=true',...
    'every y tick label/.append style={text width=width("$-0.5$"),align=right}',...
    'every y tick scale label/.style={at={(0,1)},above left,inner sep=0pt,yshift=0.3em,text width=width("$\cdot 10^{-2}$")}'...
    });

end

