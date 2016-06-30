function linkedtracedata=linkancestry_4generations(tracedata,tracestats,samplecellsID)
%links traces across 4 generations of cells. It picks the granddaughter
%with the longest trace. If you want to use different criteria to pick
%which granddauther to use, alter lines 23:29
linkedtracedata=ones(size(tracedata))*NaN;
for i=1:numel(samplecellsID)
    orgcellid=samplecellsID(i);
    cellid=orgcellid;
    condition=true;
    while condition
        goodframes=find(~isnan(tracedata(cellid,:,1)));
        linkedtracedata(orgcellid,goodframes,:)=tracedata(cellid,goodframes,:);
        cellid=tracestats(cellid,4);
        condition=~isnan(cellid);
    end
    %% add granddaughter data if it exists
    granddaughter_cells=[];
    granddaughter_cells=find(tracestats(:,4)==orgcellid);
    if ~isempty(granddaughter_cells)
        if length(granddaughter_cells)>=2
            goodframes_granddaughter_one=find(~isnan(tracedata(granddaughter_cells(1),:,1)));
            goodframes_granddaughter_two=find(~isnan(tracedata(granddaughter_cells(2),:,1)));
            if length(goodframes_granddaughter_one)>=length(goodframes_granddaughter_two)
                goodframes_granddaughter=goodframes_granddaughter_one;
                granddaughterid=granddaughter_cells(1);
            else
                goodframes_granddaughter=goodframes_granddaughter_two;
                granddaughterid=granddaughter_cells(2);
            end
        elseif length(granddaughter_cells)==1;
            goodframes_granddaughter=find(~isnan(tracedata(granddaughter_cells,:,1)));
            granddaughterid=granddaughter_cells;
        end
        linkedtracedata(orgcellid,goodframes_granddaughter,:)=tracedata(granddaughterid,goodframes_granddaughter,:);
        %% add great granddaughter data if it exists
        greatgranddaughter_cells=[];
        greatgranddaughter_cells=find(tracestats(:,4)==granddaughterid);
        if ~isempty(greatgranddaughter_cells)
            if length(greatgranddaughter_cells)>=2
                goodframes_greatgranddaughter_one=find(~isnan(tracedata(greatgranddaughter_cells(1),:,1)));
                goodframes_greatgranddaughter_two=find(~isnan(tracedata(greatgranddaughter_cells(2),:,1)));
                if length(goodframes_greatgranddaughter_one)>=length(goodframes_greatgranddaughter_two)
                    goodframes_greatgranddaughter=goodframes_greatgranddaughter_one;
                    greatgranddaughterid=greatgranddaughter_cells(1);
                else
                    goodframes_greatgranddaughter=goodframes_greatgranddaughter_two;
                    greatgranddaughterid=greatgranddaughter_cells(2);
                end
            elseif length(greatgranddaughter_cells)==1
                goodframes_greatgranddaughter=find(~isnan(tracedata(greatgranddaughter_cells,:,1)));
                greatgranddaughterid=greatgranddaughter_cells;
            end
            linkedtracedata(orgcellid,goodframes_greatgranddaughter,:)=tracedata(greatgranddaughterid,goodframes_greatgranddaughter,:);
        end
    end
end
