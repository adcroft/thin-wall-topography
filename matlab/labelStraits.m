function [] = labelstraits

[x,y,angle,label]=textread('Straits.txt','%f%f%f%s','delimiter','\t','commentstyle','matlab');

for n=1:length(x)
 if y(n)>min(ylim) && y(n)<max(ylim) && x(n)>min(xlim) && x(n)<max(xlim) 
  h=text(x(n),y(n),label{n});
  set(h,'Rotation',angle(n),'HorizontalAlignment','center')
  set(h,'FontSize',10)
  set(h,'Color','w')
 end
end

tracks={'denmark' 'faroe'};
for j=1:length(tracks)
 load(['track_' tracks{j} '.mat'])
 line(track.x,track.y,'Color','r')
end
