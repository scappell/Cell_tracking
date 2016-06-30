function tracestats = getstats(tracedata,genealogy)
numcells=size(tracedata,1);
tracestats=ones(numcells,4)*NaN;
for c=1:numcells
    tracestats(c,1)=find(~isnan(tracedata(c,:,1)),1,'first');
    tracestats(c,2)=find(~isnan(tracedata(c,:,1)),1,'last');
end
tracestats(:,3)=tracestats(:,2)-tracestats(:,1)+1;
tracestats(:,4)=genealogy;
