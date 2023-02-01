function [Dst,Dsp]=projDist(t1,t2,d1,d2,pdata)

% Dst= adjusted distance 
    % = at each time point, averaged euclidean dist of each point to all 
    % points in the other condition group, devided by the dispersion
    % measured at that point
    
% Dsp = dispersion 
    % = self distance
    % = each point's distance to all other points in the same condtion
    % group

% project
for t1n=1:length(t1)
    for t2n=1:length(t2)
        trldata=[centering(d1(:,:,t1n),1);centering(d2(:,:,t2n),1)];
        if t1n==1 && t2n==1
            top70PCs=trldata*pdata.loadings(:,1:pdata.top70);
            tmp1(:,:,1)=top70PCs(1:100,:);
            tmp2(:,:,1)=top70PCs(101:200,:);
        else
            top70PCs=trldata*pdata.loadings(:,1:pdata.top70);
            tmp1(:,:,end+1)=top70PCs(1:100,:);
            tmp2(:,:,end+1)=top70PCs(101:200,:);
        end
    end
end

% keyboard;
% pcs1=permute(tmp1,[3 2 1]);
% pcs2=permute(tmp2,[3 2 1]);
% reshape: take the first dim first and fill in the first dim first and
% then take/fill second dim, & so on so forth

tt1=reshape(permute(tmp1,[3 2 1]),[length(t1) length(t2) size(tmp1,2) size(tmp1,1)]);
tt2=reshape(permute(tmp2,[3 2 1]),[length(t1) length(t2) size(tmp2,2) size(tmp2,1)]);

pcs1=reshape(nanmean(tt1,2),[length(t1) size(tmp1,2) size(tmp1,1)]);
pcs2=reshape(nanmean(tt2,1),[length(t2) size(tmp1,2) size(tmp1,1)]);

% calculate distance
% do not plot -- too messy
for tp=1:100
    clear var1 var2
    var1=pcs1(:,:,tp);
    var2=pcs2(:,:,tp);
    mD=adjDist(var1,var2);
    Dst.aD12(:,tp)=mD.aD12;
    Dst.aD21(:,tp)=mD.aD21;
    
    Dsp.cnd1(:,tp)=EucDisperse(var1);
    Dsp.cnd2(:,tp)=EucDisperse(var2);
    
end

end