function [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_1(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite)
%%% determine median cell size for blur detection %%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=zeros(3,1); numcells=zeros(3,1);
for i=1:3
    raw1=log(single(imread([rawdir,name1,num2str(frames(i)),'.tif'])));
    nuc_mask=blobdetector_4(log(raw1),nucr,blobthreshold,debrisarea);
    [~,numcells(i)]=bwlabel(nuc_mask);
    nuc_area(i)=median(cell2mat(struct2cell(regionprops(nuc_mask,'Area'))));
end
dims=size(nuc_mask);
height=dims(1); width=dims(2);
blurthreshhigh=1.10*nanmedian(nuc_area);
%%% determine first good frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstgoodindex=find(numcells>0 & nuc_area<blurthreshhigh,1,'first');
if firstgoodindex>1
    for i=1:firstgoodindex-1
        badframes(frames(1)-1+i)=1;
        if maskwrite
            imwrite(uint16(zeros(dims)),[maskdir,namenucedge,num2str(frames(i)),'.tif']);
        end
    end
end
badframes(frames(1)-1+firstgoodindex)=0;
blurthreshhigh=1.08*nuc_area(firstgoodindex);
blurthreshlow=0.95*nuc_area(firstgoodindex);
numthresh=0.8*numcells(firstgoodindex); %was 0.5
end