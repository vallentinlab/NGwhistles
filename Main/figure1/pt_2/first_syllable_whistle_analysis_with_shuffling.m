%% Whistle analysis for PAIRS matches TAKING JUST THE FIRST SYLLABLE %%
clear all
close all
clc

[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*.xls');
%%

% for
bird=1
    

%%
filename= char(strcat(pathname,'\',filelist(bird,1).name));
[data,text] = xlsread(filename);

stim=data(:,2);
resp=data(:,3);

q=[0,9000];
qq=[0,9000];

%%
figure(bird)
subplot(2,2,1)
scatter(stim,resp)
xlim([1000 9000]);
ylim([1000 9000]);
axis square
hold on
plot (q,qq, '-k')
set(gca,'TickDir','out')
set(gca, 'xScale', 'log')
set(gca, 'yScale', 'log')

ylabel('First whistle syllable pitch of matching bird (Hz)')

xlabel('Whistle syllable pitch of matched bird (Hz)')

%% shuffled data

fake_dist=[];
fake_median=[];

fake_perc_inside_100=[];

for i=1:1000
    
fake_stim=stim(randperm(length(stim)));
fake_resp=resp(randperm(length(resp)));

XX=[fake_stim fake_resp];
d=point_to_line_distance(XX,[1000 1000 ],[ 9000 9000] );

fake_dist=[fake_dist d];
fake_median=[fake_median median(d)];



g=fake_stim-fake_resp;
[i,ii]=find(g<0);
Q=[fake_stim(i) fake_resp(i)];
D=point_to_line_distance(Q,[1000 1000 ],[ 9000 9000] );

g=fake_stim-fake_resp;
[i,ii]=find(g>0);
Q=[fake_stim(i) fake_resp(i)];
DD=point_to_line_distance(Q,[1000 1000 ],[ 9000 9000] );

DDD=[-D; DD];


[k,kk]=find(DDD<100 & DDD>-100);
%[k,kk]=find(DDD<mean(DDD)+std(DDD) & DDD>-(mean(DDD)+std(DDD)));


%[k,kk]=find(DDD<(mean(DDD)+std(DDD)/sqrt(length(DDD))) & DDD>(mean(DDD)-(std(DDD)/sqrt(length(DDD)))));

fake_perc_inside_100=[fake_perc_inside_100, 100-sum(kk)/length(g)*100];

end
figure (2)
subplot(1,3,2)
histogram(DDD, 20, 'Normalization','probability')
hold on
%xline(median(DDD), 'r-')
%xline(-(mean(DDD)+std(DDD)), 'k--')
%xline((mean(DDD)+std(DDD)), 'k--')

xline(-100, 'k-')
xline(100, 'k-')
%xline(0,'k-')
xlim([-5000 5000]);
axis square
hold on
set(gca,'TickDir','out')
box off
ylim([0 0.3]);

ylabel('Probability')

xlabel('Euclidian distance (Hz)')

title('Expected (example)')


%% real data distance

g=stim-resp;
[i,ii]=find(g<0);
Q=[stim(i) resp(i)];
D=point_to_line_distance(Q,[1000 1000 ],[ 9000 9000] );

g=stim-resp;
[i,ii]=find(g>0);
Q=[stim(i) resp(i)];
DD=point_to_line_distance(Q,[1000 1000 ],[ 9000 9000] );

DDD=[-D; DD];

[k,kk]=find(DDD<100 & DDD>-100);
%[k,kk]=find(DDD<(mean(DDD)+std(DDD)/sqrt(length(DDD))) & DDD>(mean(DDD)-(std(DDD)/sqrt(length(DDD)))));

100-sum(kk)/length(g)*100

figure (2)
subplot(1,3,3)

histogram(fake_perc_inside_100,20, 'Normalization','probability')
xline(100-sum(kk)/length(g)*100, 'r-')
%xlim([0 30 ]);
axis square
hold on
set(gca,'TickDir','out')
box off
ylabel('Probability')
title('Observed vs. Expected (1000 times)')

xlabel('Data with |Euclidian distance| larger than 100 Hz (%)')
ylim([0 0.3]);
xlim([75 100]);


subplot(1,3,1)
histogram(DDD, 20, 'Normalization','probability')
hold on
%xline(median(DDD), 'r--')
%xline(-(mean(DDD)+std(DDD)), 'k--')
%xline((mean(DDD)+std(DDD)), 'k--')

xline(-100, 'k-')
xline(100, 'k-')
%xline(0,'k-')
xlim([-5000 5000]);
axis square
hold on
set(gca,'TickDir','out')
box off
ylim([0 0.3]);
title('Observed')

ylabel('Probability')

xlabel('Euclidian distance (Hz)')


%%

figure(bird)

subplot(2,2,2)
histogram(DDD, 20, 'Normalization','probability')
hold on
xline(median(DDD), 'r-')
%xline(-(mean(DDD)+std(DDD)), 'k--')
%xline((mean(DDD)+std(DDD)), 'k--')

xline((median(DDD)-100), 'k--')
xline(median(DDD)+100, 'k--')
%xline(0,'k-')
xlim([-5000 5000]);
axis square
hold on
set(gca,'TickDir','out')
box off
ylim([0 0.3]);

ylabel('Probability')

xlabel('Euclidian distance (Hz)')

%%

X=[stim resp];
d=point_to_line_distance(X,[1000 1000 ],[ 9000 9000] );
real_median=median(d)

figure(bird)

subplot(2,2,3)
histogram(d,20,'Normalization','probability')
axis square
hold on
set(gca,'TickDir','out')
box off
xlim([0 5000]);
hold on
xline (real_median, 'r')
ylim([0 0.3]);

ylabel('Probability')

xlabel('|Euclidian distance| (Hz)')


%% figure

figure(bird)
subplot(2,2,4)
histogram(fake_median,10,'Normalization','probability')

hold on
xline (real_median, 'r')

ylim([0 0.3]);

xlim([400 1400]);
axis square
hold on
set(gca,'TickDir','out')
box off

ylabel('Probability')

xlabel('|Euclidian distance| (Hz)')


% d=stim-resp;
% figure
% histogram(d,25,'Normalization','probability')
% xline(mean(d), 'r--')
% xline(std(d), 'k--')
% xline(-std(d), 'k--')
% xlim([-7000 7000]);
% axis square
% set(gca,'TickDir','out')
% box off

% scatterDiagHist(stim, resp)
% axis square
% xlim([1000 9000]);
% ylim([1000 9000]);

% figure
% histogram(d,20,'Normalization','probability')
% hold on
% xline (real_median, 'r')
%  xlim([0 5000]);
%     axis square
% hold on
% set(gca,'TickDir','out')
% box off


%end