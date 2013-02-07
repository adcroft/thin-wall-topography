function [] = drawColorbar(Hrng,ticks)
% drawColorbar(range,ticks)
%
% Draws a stand alone colorbar with data range "range" and tick marks at "ticks"
% e.g.
% >> drawColorbar([-6500 1500],-7000:1000:7000)

% Copyright 2008-2013 Alistair Adcroft, Princeton University.
%
% This file is part of the thin-wall-topography software suite.
%
% thin-wall-topography is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% thin-wall-topography is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with Foobar. If not, see http://www.gnu.org/licenses/.

pwidth=4;
set(gcf,'PaperSize',[8.5 1]*pwidth/8.5,'PaperPosition',[0 0 8.5 1]*pwidth/8.5)
ddH=diff(Hrng)/length(colormap)/10;
HHx=(min(Hrng)-diff(Hrng)/400/2):50:max(Hrng);
HHx=(0:10*length(colormap))*ddH+min(Hrng);
pcol(HHx,[0 1],[HHx;HHx]); caxis(Hrng);
set(gca,'Position',[0.02 0.35 0.96 0.61],'YTick',[])
for l=length(ticks)
% if tick0s(l)>min(Hrng) && ticks(l)<max(Hrng)
    line([1 1]*ticks(l),[-1 2],'Color','k');
% end
end
%set(gca,'XTick',ticks)
