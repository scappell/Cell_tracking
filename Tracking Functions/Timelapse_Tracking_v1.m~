function Timelapse_Tracking_v1(row,col,site,settings,debug_mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if debug_mode
    settings.EndFrame=settings.StartFrame+2;
end
if debug_mode
    keyboard;
end
%%% Designate which Well is being analyzied %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datadir=settings.projectpath;
imagepath=settings.imagepath;
shadingpath=settings.shadingpath;
if ~exist(datadir,'dir')
    mkdir(datadir);
end
if ~exist([datadir,'tracedata_',shot,'.mat']);
    separatedirectories=settings.separatedirectories_option;
    if separatedirectories==1
        rawdir=[imagepath,'Raw/',shot,'/',shot,'_'];
        maskdir=[imagepath,'Mask/',shot];
        biasdir=[imagepath,'Bias/'];
    else
        rawdir=[imagepath,'Raw/',shot,'_'];
        maskdir=[imagepath,'Mask'];
        biasdir=[imagepath,'Bias/'];
    end
    maskwrite=settings.maskwrite_option;
    if ~exist(maskdir,'dir') && maskwrite
        mkdir(maskdir);
    end
    if separatedirectories==0
        maskdir=[maskdir,'/',shot,'_'];
    end
    %%% General settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SF=settings.StartFrame;EF=settings.EndFrame;
    name1=settings.nucleus_name; %nuclear channel
    name2=settings.signal2;
    if settings.signal3_option
        name3=settings.signal3;
    end
    %%% Segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucr=settings.nucr;
    debrisarea=settings.debrisarea;
    boulderarea=settings.boulderarea;
    blobthreshold=settings.blobthreshold;
    %%% Tracking parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxjump=settings.maxjump;
    masschangethreshold=settings.masschangethreshold;
    areachangethreshold=settings.areachangethreshold;
    daughtervariance=settings.daughtervariance;
    trackparams={nucr,maxjump,debrisarea,masschangethreshold,areachangethreshold,daughtervariance};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frames=SF:EF;
    totalframes=numel(frames);
    badframes=ones(EF,1)*NaN;
    if SF>1
        badframes(1:SF-1)=0;
    end
    jitters=zeros(EF,2);
    blocksize=settings.blocksize;
    maxcellnum=blocksize;
    parameternum=11;
    tracedata=ones(maxcellnum,EF,parameternum)*NaN;
    tracking=ones(maxcellnum,5)*NaN;
    timetotal=tic;
     [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_4(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
    regheight=1:0.5*height; regwidth=1:0.5*width;
    %%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([shadingpath,'BG_bin1','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
    if settings.bias_sig1_option
        load([biasdir,name1,num2str(site),'.mat']); bias1=bias;
    end
    if settings.bias_sig2_option
        load([biasdir,name2,num2str(site),'.mat']); bias2=bias;
    end
    if settings.bias_sig3_option & settings.signal3_option
        load([biasdir,name3,num2str(site),'.mat']); bias3=bias;
    end
    
    %%% Imaging processing, Frame by Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=firstgoodindex:totalframes
        f=frames(i); fprintf('frame %0.0f\n',f);
        timeframe=tic;
        %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        raw1=double(imread([rawdir,name1,num2str(f),'.tif'])); if settings.shadingcorrection; raw1=(raw1-bgcmos);end; if settings.bias_sig1_option; raw1=raw1./bias1;end;
        raw2=double(imread([rawdir,name2,num2str(f),'.tif'])); if settings.shadingcorrection;raw2=(raw2-bgcmos);end; if settings.bias_sig2_option; raw2=raw2./bias2;end;
        if settings.signal3_option
            raw3=double(imread([rawdir,name3,num2str(f),'.tif'])); if settings.shadingcorrection;raw3=(raw3-bgcmos);end; if settings.bias_sig3_option; raw3=raw3./bias3;end;
        end
        %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i==firstgoodindex
            firstsegmethod=settings.firstsegmethod;
            switch firstsegmethod
                case 'log'
                    nuc_mask=blobdetector_4(log(raw1),nucr,blobthreshold,debrisarea);
                case 'single'
                    blurradius=settings.blurradius;
                    nuc_mask=threshmask(raw1,blurradius);
                    nuc_mask=markershed(nuc_mask,round(nucr*2/3));
                case 'double'
                    blurradius=settings.blurradius;
                    nuc_mask=threshmask(raw1,blurradius);
                    nuc_mask=markershed(nuc_mask,round(nucr*2/3));
                    nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea*2);
            end
            foreground=nuc_mask;
            nuc_mask=bwareaopen(nuc_mask,debrisarea);
            %%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
            eccentricitythresh=settings.eccentricitythresh;
            nuc_mask=excludelargeandwarped_3(nuc_mask,boulderarea,eccentricitythresh);
        else
            nuc_mask=threshmask(raw1,1);
            diagnostic_temp=bwareaopen(nuc_mask,debrisarea);
            diagnostic_cc=bwconncomp(diagnostic_temp);
            diagnostic_data=regionprops(diagnostic_cc,'basic');
            for numobjects_diagnostic=1:size(diagnostic_data,1)
                object_area(numobjects_diagnostic)=diagnostic_data(numobjects_diagnostic).Area;
            end
            if diagnostic_cc.NumObjects<=5 | max(object_area)>(size(raw1,1)*size(raw1,2))/2;
                %%% patch for bubble
                perc20=round(prctile(raw1(:),20));
                raw1(raw1<perc20)=perc20;
                nuc_mask=threshmask(raw1,1);
            end
            %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lastgoodframe=find(badframes==0,1,'last');
            [reljitx,reljity]=registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
            jitters(f,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
            secondsegmethod=settings.secondsegmethod;
            switch secondsegmethod
                case 'log'
                    nuc_mask=blobdetector_4(log(raw1),nucr,blobthreshold,debrisarea);
                case 'double'
                    nuc_mask=threshmask(raw1,blurradius);
                    nuc_mask=markershed(nuc_mask,round(nucr*2/3));
                    nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea*2);
                case 'apriori'
                    [nuc_mask,marker_mask]=apriori_markermask(nuc_mask,nuc_center,jitters(f,:));
            end
            foreground=nuc_mask;
            nuc_mask=bwareaopen(nuc_mask,debrisarea);
        end
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask=imclearborder(nuc_mask);
%%% Only for DeBug Mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if debug_mode
            keyboard;
        end
        %{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
hist(nuc_area,100);

anti_mask=bwareaopen(nuc_mask,debrisarea);
temp_mask=nuc_mask-anti_mask;
extractmask=bwmorph(temp_mask,'remove');

anti_mask=bwareaopen(nuc_mask,1000);
extractmask=bwmorph(anti_mask,'remove');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        compression=settings.compression;
        nanmask=imdilate(foreground,strel('disk',nucr/2));
        nanmaskcyto=imdilate(foreground,strel('disk',nucr*2));
        blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
        blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
        if settings.signal3_option
            blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
        end
        real1=bgsubmasked_global_2(blur1,nanmask,11,compression,50);
        if settings.ringcalc_sig2
            real2=bgsubmasked_global_2(blur2,nanmaskcyto,1,compression,10);
        else
            real2=imtophat(blur2,strel('disk',15,0));
        end
        if settings.signal3_option
            if settings.ringcalc_sig3
                real3=bgsubmasked_global_2(blur3,nanmaskcyto,1,compression,10);
            else
                real3=imtophat(blur3,strel('disk',15,0));
            end
        end
        %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [nuc_label,numcells]=bwlabel(nuc_mask);
        nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
        nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
        %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mednuc=median(nuc_area);
        if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)
            fprintf('badframe: frame %0.0f\n',f);
            badframes(f)=1;
            badextractmask=bwmorph(nuc_mask,'remove');
            if maskwrite
                imwrite(uint16(extractmask),[maskdir,namenucedge,num2str(f),'.tif']);
            end
            continue;
        end
        blurthreshhigh=1.1*mednuc;
        blurthreshlow=0.8*mednuc;
        numthresh=0.5*numcells;
        nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
        nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
        %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mass=nuc_density.*nuc_area;
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i>firstgoodindex
            nuc_center(:,1)=nuc_center(:,1)+jitters(f,1);
            nuc_center(:,2)=nuc_center(:,2)+jitters(f,2);
            %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
            debugpackage={extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
            %%% track & correct merges (update centers, masses and labels) %%%%
            [tracedata,curdata,tracking,nuc_label]=adaptivetrack_9(f,lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,jitters(f,:),trackparams,debugpackage);
            badframes(f)=0;
        end
        %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        extractmask=bwmorph(nuc_label,'remove');
        if maskwrite
            imwrite(uint16(extractmask),[maskdir,namenucedge,num2str(f),'.tif']);
        end
        %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cellid=find(~isnan(curdata(:,1)));
        numlivecells=numel(cellid);
        curdata=curdata(cellid,:);
        nuc_center=curdata(:,[1 2]);
        nuc_area=curdata(:,3);
        nuc_mass=curdata(:,4);
        nuc_info=regionprops(nuc_label,'PixelIdxList');
        nanvec=ones(numlivecells,1)*NaN; sig1=nanvec; sig2=nanvec; sig3=nanvec;
        sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec;
        sig3ring_75th=nanvec; sig3ring_fgmedian=nanvec;
        for n=1:numlivecells
            cc=cellid(n);
            sig1(n)=mean(real1(nuc_info(cc).PixelIdxList));
            sig2(n)=median(real2(nuc_info(cc).PixelIdxList));
            if settings.signal3_option
                sig3(n)=median(real3(nuc_info(cc).PixelIdxList));
            else
                sig3(n)=median(real2(nuc_info(cc).PixelIdxList));
            end
        end
        if settings.ringcalc_sig2==1 | settings.ringcalc_sig3==1
            if (settings.magnification==10 & settings.binsize==1) | (settings.magnification==20 & settings.binsize==2)
                innerrad=1; outerrad=5;
            else
                innerrad=1; outerrad=5;
            end
            ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,real2);
            ring_info=regionprops(ring_label,'PixelIdxList');
            
            sig2thresh=settings.ringthresh;
            sig3thresh=settings.ringthresh;
            for n=1:numlivecells
                cc=cellid(n);
                if cc>numel(ring_info)
                    break;
                end
                ring2all=real2(ring_info(cc).PixelIdxList);
                ring2all(ring2all>prctile(ring2all,98))=[];
                sig2ring_75th(n)=prctile(ring2all,75);
                ring2foreground=ring2all(ring2all>sig2thresh);
                if numel(ring2foreground)<100
                    ring2foreground=ring2all;
                end
                if numel(ring2all)>100
                    sig2ring_fgmedian(n)=nanmedian(ring2foreground);
                end
                
                if settings.signal3_option
                    ring3all=real3(ring_info(cc).PixelIdxList);
                    ring3all(ring3all>prctile(ring3all,98))=[];
                    sig3ring_75th(n)=prctile(ring3all,75);
                    ring3foreground=ring3all(ring3all>sig3thresh);
                    
                    if numel(ring3foreground)<100
                        ring3foreground=ring3all;
                    end
                    if numel(ring3all)>100
                        sig3ring_fgmedian(n)=nanmedian(ring3foreground);
                    end
                end
            end
        end
        
        %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig3,sig3ring_75th,sig3ring_fgmedian];
        if maxcellnum-max(cellid)<blocksize
            tempdata=ones(blocksize,EF,parameternum)*NaN;
            temptrack=ones(blocksize,5)*NaN;
            tracedata=[tracedata;tempdata];
            tracking=[tracking;temptrack];
            maxcellnum=maxcellnum+blocksize;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toc(timeframe);
    end
        [tracedata,genealogy,jitters]=postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
        %%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
        toc(timetotal);
        clear all; clear mex;
end
clear all; clear mex;