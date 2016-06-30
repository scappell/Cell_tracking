function vIdx=getdeflections(orderedset,nucr)

%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
offsetshort=nucr/4;
if offsetshort==1
    offsetshort=2;
end
offsetlong=2*offsetshort;
gradientoffsetshort=offsetshort;
gradientoffsetlong=offsetlong;
shortgradthresh = pi/6;
longgradthresh = pi/6;

%%% calculate angular deflection at each point of the boundary %%%%%%%%%%%%
perilength=size(orderedset,1);
%%%%%%% short steps %%%%%%%%%%%%%%%%%%
orderedsetoffsetshort=[orderedset(offsetshort+1:end,:);orderedset(1:offsetshort,:)];
shortdiff=orderedsetoffsetshort-orderedset;
shortgrad=atan2(shortdiff(:,2),shortdiff(:,1));   %angle in radians
shortgradoffset=[shortgrad(gradientoffsetshort+1:end,:);shortgrad(1:gradientoffsetshort,:)];
shortgraddiff=shortgradoffset-shortgrad;
shortgraddiff=shortgraddiff+2*pi*(shortgraddiff<0);  %account for 4 quadrants
%%%%%%% long steps %%%%%%%%%%%%%%%%%%%
orderedsetoffsetlong=[orderedset(offsetlong+1:end,:);orderedset(1:offsetlong,:)];
longdiff=orderedsetoffsetlong-orderedset;
longgrad=atan2(longdiff(:,2),longdiff(:,1));
longgradoffset=[longgrad(gradientoffsetlong+1:end,:);longgrad(1:gradientoffsetlong,:)];
longgraddiff=longgradoffset-longgrad;
longgraddiff=longgraddiff+2*pi*(longgraddiff<0);

%%% find deflections above threshold %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% short steps %%%%%%%%%%%%%%%%%%
shortgraddiff(shortgraddiff>=pi)=0;  %exclude
vIdxmaskshort=shortgraddiff>shortgradthresh;
%%%%%%% long steps %%%%%%%%%%%%%%%%%%%
vIdxmasklong=longgraddiff>longgradthresh & longgraddiff<pi;
vIdxmasklong=[zeros(offsetlong,1);vIdxmasklong(1:end-offsetlong)];
vIdxmasklong=imdilate(vIdxmasklong,strel('square',1+nucr/2));

%%% find local maxima of short steps %%%%%%%%%%%%%%%%%%%%%%%
vIdxmaskshort=imclose(vIdxmaskshort,strel('square',3));  %join proximal deflection islands
vIdxobs=regionprops(bwlabel(vIdxmaskshort),'PixelIdxList');
maxmask=zeros(size(vIdxmaskshort));
for rpc=1:length(vIdxobs)
    pix=vIdxobs(rpc).PixelIdxList;
    [~,index]=max(shortgraddiff(pix));
    maxmask(pix(index)+offsetshort)=1;
end
maxmask=maxmask(1:perilength);  %remove any overhang

%%% find coincidence of long mask & local maxima of short mask %%%%%%%%%%%%
vIdxmask=vIdxmasklong & maxmask;
vIdx=find(vIdxmask);

end