function DistByTimePlot3_e2(indata,pccdata,outdata,smo,yLims) %,yLims

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

% keyboard

inMean=nanmean(indata,1);
inSem=confidenceInterval(indata,1);

plot(smooth(inMean,smo),'r-','LineWidth',1);
hold on
ylim(yLims)
plot(smooth(inMean-inSem,smo),'r:','LineWidth',1);
plot(smooth(inMean+inSem,smo),'r:','LineWidth',1);

pccMean=nanmean(pccdata,1);
pccSem=confidenceInterval(pccdata,1);

plot(smooth(pccMean,smo),'b-','LineWidth',1);
plot(smooth(pccMean-pccSem,smo),'b:','LineWidth',1);
plot(smooth(pccMean+pccSem,smo),'b:','LineWidth',1);


outMean=nanmean(outdata,1);
outSem=confidenceInterval(outdata,1);

plot(smooth(outMean,smo),'k-','LineWidth',1);
plot(smooth(outMean-outSem,smo),'k:','LineWidth',1);
plot(smooth(outMean+outSem,smo),'k:','LineWidth',1);

vline(11);
vline(31);
vline(41);
vline(51);
vline(71);
vline(91);
legend('In','','','PCC','','','Out','','',...
    num2str(binStarts(11)),num2str(binStarts(31)),...
    num2str(binStarts(41)),num2str(binStarts(51)),...
    num2str(binStarts(71)),num2str(binStarts(91)));
xlabel('Time');
ylabel('Distance');


end