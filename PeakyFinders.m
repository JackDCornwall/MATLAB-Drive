close all;
clear all;

%%%% Read in the data %%%%

Data=dlmread('test_peaky.txt');

t = Data(:,1)/1000;
yin = Data(:,2);

%%%% Plot the oscillation data %%%%

figure;
scatter(t,yin,'filled');
hold;
plot(t,yin);
set(gca,'FontSize',14);
xlabel('time (s)');
ylabel('Pendulum angle (arb. units)');

%%
%%%% Finding peaks in the oscillation data %%%%
ndata = length(yin);
offset = mean(yin((ndata-100):ndata));

peakcount=0;
for i=5:length(t)-5
    grad1=(yin(i)-yin(i-4))/(t(i)-t(i-4));
    grad2=(yin(i+4)-yin(i))/(t(i+4)-t(i));
    if grad1 > 0 && grad2 < 0
        peakcount=peakcount+1;
        peakpos(peakcount)=t(i);
        peakheight(peakcount)=yin(i)-offset;
    end
end

% %%%% Plot the peak positions vs peak heights %%%%
% figure;
% scatter(peakpos,peakheight,'filled');
% hold;

npeaks = int64(0.5*length(peakpos));

%%%% Chop out the first half of the peaks %%%%
peakpos=peakpos(1:npeaks);
peakheight=peakheight(1:npeaks);
PeakHeightError=ones(1,npeaks);
PeakHeightError=PeakHeightError/2;
logpeakheight=log(abs(peakheight));
logpeakheighterror=abs(PeakHeightError./peakheight);

OutData = [peakpos' logpeakheight' logpeakheighterror'];

dlmwrite('OutputDecayData.txt',OutData,'\t');

figure;
scatter(t,yin,'filled');
hold;
scatter(peakpos,peakheight+offset,'filled','r');
set(gca,'FontSize',14);
xlabel('time (s)');
ylabel('Pendulum angle (arb. units)');