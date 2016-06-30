function avgtheta=getslope_forward_avg(series,step)
numsteps=length(step);
thetalength=length(series)-max(step);
theta=zeros(numsteps,thetalength);
for cc=1:numsteps
    i=step(cc);
    temp=atan2(series(1+i:end)-series(1:end-i),i);
    theta(cc,:)=temp(1:thetalength);
end
%avgtheta=sum(theta,1);
avgtheta=mean(theta,1);
avgtheta=[avgtheta,ones(1,max(step))*avgtheta(end)];