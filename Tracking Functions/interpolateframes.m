function [tracedata,jitters]=interpolateframes(tracedata,jitters,badframes,tracking)
totalframes=size(tracedata,2);
badframeid=find(badframes);
for i=badframeid'
    if i==1
        firstgood=find(badframes==0,1,'first');
        tracedata(:,i,:)=tracedata(:,firstgood,:);
        jitters(i,:)=jitters(firstgood,:);
    elseif i==totalframes || badframes(i+1)==1
        tracedata(:,i,:)=tracedata(:,i-1,:);
        jitters(i,:)=jitters(i-1,:);
    else
        prevnum=find(~isnan(tracedata(:,i-1,1)),1,'last');
        tracedata(1:prevnum,i,:)=(tracedata(1:prevnum,i-1,:)+tracedata(1:prevnum,i+1,:))/2;
        jitters(i,:)=mean([jitters(i-1,:);jitters(i+1,:)]);
    end
    badframes(i)=0;
end

daughters=find(~isnan(tracking(:,1)))';
for d=daughters
    m=tracking(d,1);
    daughterfirst=find(~isnan(tracedata(d,:,1)),1,'first');
    motherlast=find(~isnan(tracedata(m,:,1)),1,'last');
    if daughterfirst==motherlast+2 %badframe occurred at mitosis
        tracedata(d,daughterfirst-1,:)=(tracedata(m,motherlast,:)+tracedata(d,daughterfirst,:))/2;
    end
end

end