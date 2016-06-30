function mask=segmentdeflections_bwboundaries(mask,nucr,debrisarea)
[B,L]=bwboundaries(mask,'noholes');
obnum=max(L(:));
bordermask=zeros(size(mask));
for ci=1:obnum
    orderedset=B{ci};
    % If using bwboundaries, reverse the order.
    orderedset=[orderedset(end:-1:1,2) orderedset(end:-1:1,1)];
    bordermask=splitdeflections_4_bwboundaries(orderedset,bordermask,nucr);
end
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
mask=bwareaopen(mask,debrisarea);
end