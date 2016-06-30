function [new_mask,marker_mask]=apriori_markermask(cur_mask,prev_center,jitters)
[height,width]=size(cur_mask);
markerx=round(prev_center(:,1)-jitters(1));
markery=round(prev_center(:,2)-jitters(2));
markernan=isnan(markerx) | markerx<1 | markerx>width | markery<1 | markery>height;
markerx(markernan)=[]; markery(markernan)=[];
markeridx=sub2ind([height,width],markery,markerx);
marker_mask=zeros(height,width);
marker_mask(markeridx)=1; marker_mask=imdilate(marker_mask,strel('disk',4,0));
[marker_label,markernum]=bwlabel(marker_mask);
maskedmarkerlabels=marker_label.*cur_mask;
uniquemarkerlabels=unique(maskedmarkerlabels); uniquemarkerlabels(1)=[];
missinglabels=find(~ismember(1:markernum,uniquemarkerlabels));
marker_mask(ismember(marker_label,missinglabels))=0;
new_mask=markershed_apriori(cur_mask,marker_mask);
end