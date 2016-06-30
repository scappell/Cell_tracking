function blockimg=blockpercentile_blockimage(initimg,blocknum,prctilethresh)
getpercentile=@(block_struct) prctile(block_struct.data(:),prctilethresh);
[height,width]=size(initimg);
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
blockimg=blockproc(initimg,[blockheight blockwidth],getpercentile);
%prctileimage=imresize(blockimg,[height width],'bicubic');
end