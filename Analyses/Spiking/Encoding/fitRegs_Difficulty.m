%%
function md = fitRegs_Difficulty(data)
% for figure 3A-F: t = 435; binSize = 25;
cutSafe = 1;
for i=1:length(data)
    
    c=data(i);
    
    if cutSafe
        temp=find(c.vars(:,9)~=0);% 0:Safe 1:Lose 2:Win
        psth=c.psth(temp,:);
        vars=c.vars(temp,:);
    else
        psth=c.psth;
        vars=c.vars;
    end
    
    % NORMALIZE
    psth=psth-mean(psth(:));
    psth=psth./var(psth(:));
    
    lp=vars(:,3);% prob of left opt
    temp=vars(:,6);% magnitude of left opt
    templ=temp+1;
    rwdsize=[165 240 125];
    lm=rwdsize(templ)';
    clear temp;
    
    rp=vars(:,4);% prob of right opt
    temp=vars(:,7);% magnitude of right opt
    tempr=temp+1;
    rm=rwdsize(tempr)';
    clear temp;
    
    firstAppear=vars(:,5);% L=1; R=2
    choice=vars(:,8); % 1:Left 2:Right
    choLR=choice; % 1:Left 2:Right
    choLR(choLR==2)=-1;% 1=left -1=right
    outcome=vars(:,9);% 0:Safe 1:Lose 2:Win
    outcome(outcome==1)=0;
    outcome(outcome==2)=1;
    
    lFirstTrl=find(firstAppear==1);
    rFirstTrl=find(firstAppear==2);
    
    
    p1=zeros([length(vars) 1]);
    p2=zeros([length(vars) 1]);
    m1=zeros([length(vars) 1]);
    m2=zeros([length(vars) 1]);
    
    
    p1(lFirstTrl)=lp(lFirstTrl);
    p1(rFirstTrl)=rp(rFirstTrl);
    % p1=normalize(p1,'center','mean');
    
    
    p2(lFirstTrl)=rp(lFirstTrl);
    p2(rFirstTrl)=lp(rFirstTrl);
    
    m1(lFirstTrl)=lm(lFirstTrl);
    m1(rFirstTrl)=rm(rFirstTrl);
    
    m2(lFirstTrl)=rm(lFirstTrl);
    m2(rFirstTrl)=lm(rFirstTrl);
    
    ev1=p1.*m1;
    ev2=p2.*m2;
    
    choOpt=[];
    choOpt(1:length(ev1),1)=2;
    choOpt(choice==firstAppear)=1;
    
    winTrl=find(outcome~=1);
    cho1=find(firstAppear==choice);
    rej1=find(firstAppear~=choice);
    
    % checks
    %     mean(firstAppear==2)
    %     length(cho1Trl)/length(firstAppear)
    %     intersect(cho1,rej1) %should be empty
    
    lev=lp.*lm;
    rev=rp.*rm;
    lCorr=intersect(find(lev>rev),find(choice==1));
    rCorr=intersect(find(lev<rev),find(choice==2));
    correctTrls=sort([lCorr; rCorr],'ascend');
    numCorrect(i)=length(correctTrls);
    totalTrl(i)=length(vars(:,1));
    
    receivedRwd=outcome;
    receivedRwd(choOpt==1)=m1(choOpt==1).*outcome(choOpt==1);
    receivedRwd(choOpt==2)=m2(choOpt==2).*outcome(choOpt==2);
    
    cho1Trl=intersect(cho1,winTrl);
    rej1Trl=intersect(rej1,winTrl);
    
    %   normalize(A,2,'center','mean')
    %   normalize each row across the 2nd dim col
    
    m1=zscore(m1);
    m2=zscore(m2);
    p1=zscore(p1);
    p2=zscore(p2);
    ev1=zscore(ev1);
    ev2=zscore(ev2);
    lev=zscore(lev);
    rev=zscore(rev);
    rwdOtc=zscore(receivedRwd);
    
%     keyboard;
    difficulty = abs(ev1-ev2);
    cutoff = median(difficulty);
    difficult = ones([length(ev1),1]);
    difficult(difficulty>cutoff) = 0;
    difficult = logical(difficult);
    
    choOpt=categorical(choOpt);
    choLR=categorical(choLR);
    firstAppear=categorical(firstAppear);
    
    % atOFFER1: opt1=301:400;opt2=401:500; chocie: 500:614; outcome: 615:730 (according to trial length)
    
    % use only 400:700 [offer2 + choice]

    
    slideSize=20; % 200 ms
    useBins = 400:700;
    for t=1:length(useBins)
        
        subpsth=psth(:,useBins(t):useBins(t)+slideSize);
        slideFR=mean(subpsth,2);
        
        % EV1 vs. EV2
        md16=fitlm([ev1(difficult), ev2(difficult)], slideFR(difficult));
        
        md(1).hard.bev1(i,t)=md16.Coefficients.Estimate(2);
        md(1).hard.tev1(i,t)=md16.Coefficients.tStat(2);
        md(1).hard.pev1(i,t)=md16.Coefficients.pValue(2);
        
        md(1).hard.bev2(i,t)=md16.Coefficients.Estimate(3);
        md(1).hard.tev2(i,t)=md16.Coefficients.tStat(3);
        md(1).hard.pev2(i,t)=md16.Coefficients.pValue(3);
        
        clear md16
        md16=fitlm([ev1(~difficult), ev2(~difficult)], slideFR(~difficult));
        
        md(1).easy.bev1(i,t)=md16.Coefficients.Estimate(2);
        md(1).easy.tev1(i,t)=md16.Coefficients.tStat(2);
        md(1).easy.pev1(i,t)=md16.Coefficients.pValue(2);
        
        md(1).easy.bev2(i,t)=md16.Coefficients.Estimate(3);
        md(1).easy.tev2(i,t)=md16.Coefficients.tStat(3);
        md(1).easy.pev2(i,t)=md16.Coefficients.pValue(3);
        
        
        
        % full model 3
        tbl=table(lev(difficult),rev(difficult),choOpt(difficult),choLR(difficult),rwdOtc(difficult),firstAppear(difficult),slideFR(difficult),...
            'VariableNames',{'lev','rev','choOpt','choLR','rwdOtc','firstAppear','slideFR'});
        md17=fitlm(tbl, 'slideFR~lev+rev+choOpt+choLR+rwdOtc+firstAppear');
        
       
        md(2).hard.blev(i,t)=md17.Coefficients.Estimate(2);
        md(2).hard.tlev(i,t)=md17.Coefficients.tStat(2);
        md(2).hard.plev(i,t)=md17.Coefficients.pValue(2);
        
        md(2).hard.brev(i,t)=md17.Coefficients.Estimate(3);
        md(2).hard.trev(i,t)=md17.Coefficients.tStat(3);
        md(2).hard.prev(i,t)=md17.Coefficients.pValue(3);
        
        
        % full model 3
        clear tbl md17
        tbl=table(lev(~difficult),rev(~difficult),choOpt(~difficult),choLR(~difficult),rwdOtc(~difficult),firstAppear(~difficult),slideFR(~difficult),...
            'VariableNames',{'lev','rev','choOpt','choLR','rwdOtc','firstAppear','slideFR'});
        md17=fitlm(tbl, 'slideFR~lev+rev+choOpt+choLR+rwdOtc+firstAppear');
        
        
        
        md(2).easy.blev(i,t)=md17.Coefficients.Estimate(2);
        md(2).easy.tlev(i,t)=md17.Coefficients.tStat(2);
        md(2).easy.plev(i,t)=md17.Coefficients.pValue(2);
        
        md(2).easy.brev(i,t)=md17.Coefficients.Estimate(3);
        md(2).easy.trev(i,t)=md17.Coefficients.tStat(3);
        md(2).easy.prev(i,t)=md17.Coefficients.pValue(3);
        
        
        clear tbl
        
        disp(['Cell  ' num2str(i) '  Slide  ' num2str(t)])
    end
end

end


