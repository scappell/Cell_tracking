function [bx,by] = bridge(vpos1,vpos2)
lengthx=vpos2(1)-vpos1(1);
lengthy=vpos2(2)-vpos1(2);
longerside=max([abs(lengthx) abs(lengthy)]);
stepx=lengthx/longerside;
stepy=lengthy/longerside;
bx=zeros(longerside+1,1);by=zeros(longerside+1,1);
for bs=0:longerside
    bx(1+bs)=vpos1(1)+round(bs*stepx);
    by(1+bs)=vpos1(2)+round(bs*stepy);
end
end