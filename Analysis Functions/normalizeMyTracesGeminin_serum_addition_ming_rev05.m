function gemininNorm=normalizeMyTracesGeminin_serum_addition_ming_rev05(geminin,basal1)
%Normalizes traces by dividing by maximum

%default input
if nargin<2
    basal1=0.01;
end

[n,p]=size(geminin);
geminin_basal=repmat(basal1,n,p);
% geminin_basal_corrected=geminin+geminin_basal;
geminin_max=max(geminin,[],2);
geminin_max_pop=prctile(geminin_max,98);
geminin_max_cutoff=prctile(geminin_max,98);
for i=1:size(geminin,1);
   if geminin_max(i)<geminin_max_cutoff;
       geminin_max(i)=geminin_max_pop;
   end
end
geminin_max_repmat=repmat(geminin_max,1,p);
gemininNorm1=geminin./geminin_max_repmat;
gemininNorm=gemininNorm1+geminin_basal;
end

