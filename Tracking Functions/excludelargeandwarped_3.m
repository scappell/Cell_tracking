function nuc_mask=excludelargeandwarped(nuc_mask,boulderarea,soliditythresh)
antimask=bwareaopen(nuc_mask,boulderarea);
nuc_mask=nuc_mask-antimask;
nuc_label=bwlabel(nuc_mask);
nuc_solidity=cell2mat(struct2cell(regionprops(nuc_label,'Solidity')));
warpedobjects=find(nuc_solidity<soliditythresh); %default 10xbin1:0.9
nuc_label(ismember(nuc_label,warpedobjects))=0;
nuc_mask=nuc_label>0;