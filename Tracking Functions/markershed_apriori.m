function maskwatershed=markershed_apriori(mask,basins)
valleys=-bwdist(~mask);
%basins=imerode(mask,strel('disk',eroderadius,0));
bigvalleys=bwdist(mask);
outerridges=watershed(bigvalleys);
outerridges=outerridges==0;
finalvalleys=imimposemin(valleys,basins | outerridges);
finalridges=watershed(finalvalleys);
mask(finalridges==0)=0;
maskwatershed=mask;
end