clear all
close all
clc

%% Load data/ adress folder %%
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);

%% choose .csv files
filelist_pairs = dir('*Songs.csv');
filelist_pre = dir('*Precontrol.csv');
filelist_pb = dir('*Playback.csv');
filelist_post = dir('*Postcontrol.csv');

%% read tables
filename_pairs = char(filelist_pairs(1,1).name);
[numbers_pairs, text_pairs] = xlsread(filename_pairs);

filename_pre = char(filelist_pre(1,1).name);
[numbers_pre, text_pre] = xlsread(filename_pre);

filename_pb = char(filelist_pb(1,1).name);
[numbers_pb, text_pb] = xlsread(filename_pb);

filename_post = char(filelist_post(1,1).name);
[numbers_post, text_post] = xlsread(filename_post);

%% percentages concatenating set 1 set 2

% w_percentage_pairs=(numbers_pairs(:,2)./numbers_pairs(:,1))*100;
% 
% w_percentage_pre1=(numbers_pre(:,2)./numbers_pre(:,1))*100;
% w_percentage_pre2=(numbers_pre(:,4)./numbers_pre(:,3))*100;
% w_percentage_pre=[w_percentage_pre1;w_percentage_pre2];
% 
% w_percentage_pb1=(numbers_pb(:,2)./numbers_pb(:,1))*100;
% w_percentage_pb2=(numbers_pb(:,4)./numbers_pb(:,3))*100;
% w_percentage_pb=[w_percentage_pb1;w_percentage_pb2];
% 
% w_percentage_post1=(numbers_post(:,2)./numbers_post(:,1))*100;
% w_percentage_post2=(numbers_post(:,4)./numbers_post(:,3))*100;
% w_percentage_post=[w_percentage_post1;w_percentage_post2];
% 
% avg_w_percentage_pairs=mean(w_percentage_pairs)
% avg_w_percentage_pre=mean(w_percentage_pre, 'omitnan')
% avg_w_percentage_pb=mean(w_percentage_pb, 'omitnan')
% avg_w_percentage_post=mean(w_percentage_post, 'omitnan')
% 
% [p,h]=ranksum(w_percentage_pairs,w_percentage_pre)
% [p,h]=ranksum(w_percentage_pre,w_percentage_pb)
% [p,h]=ranksum(w_percentage_pb,w_percentage_post)
% [p,h]=ranksum(w_percentage_pre,w_percentage_post)

%% Percentages of avg set 1 set 2

w_match_perc2=(numbers_pb(:,7)./numbers_pb(:,6))*100;
w_match_perc1=(numbers_pb(:,3)./numbers_pb(:,2))*100;
w_match_perc=mean([w_match_perc1,w_match_perc2],2,'omitnan');

 w_match_perc_over_tot1=(numbers_pb(:,3)./numbers_pb(:,1))*100;
 w_match_perc_over_tot2=(numbers_pb(:,7)./numbers_pb(:,5))*100;
w_match_perc_over_tot=mean([w_match_perc_over_tot1,w_match_perc_over_tot2],2,'omitnan');


w_percentage_pairs=(numbers_pairs(:,2)./numbers_pairs(:,1))*100;

w_percentage_pre1=(numbers_pre(:,2)./numbers_pre(:,1))*100;
w_percentage_pre2=(numbers_pre(:,4)./numbers_pre(:,3))*100;
w_percentage_pre=mean([w_percentage_pre1,w_percentage_pre2],2,'omitnan');

w_percentage_pb1=(numbers_pb(:,2)./numbers_pb(:,1))*100;
w_percentage_pb2=(numbers_pb(:,6)./numbers_pb(:,5))*100;
w_percentage_pb=mean([w_percentage_pb1,w_percentage_pb2],2,'omitnan');

w_percentage_post1=(numbers_post(:,2)./numbers_post(:,1))*100;
w_percentage_post2=(numbers_post(:,4)./numbers_post(:,3))*100;
w_percentage_post=mean([w_percentage_post1,w_percentage_post2],2,'omitnan');

avg_w_percentage_pairs=mean(w_percentage_pairs);
avg_w_percentage_pre=mean(w_percentage_pre, 'omitnan');
avg_w_percentage_pb=mean(w_percentage_pb, 'omitnan');
avg_w_percentage_post=mean(w_percentage_post, 'omitnan');

[p,h]=ranksum(w_percentage_pairs,w_percentage_pre)
[p,h]=signrank(w_percentage_pre,w_percentage_pb)
[p,h]=signrank(w_percentage_pb,w_percentage_post)
[p,h]=signrank(w_percentage_pre,w_percentage_post)
[p,h]=ranksum(w_percentage_pairs,w_percentage_pb)
[p,h]=ranksum(w_percentage_pairs,w_percentage_post)


%% Perc. of whistle songs: Pairs VS Playback experiment

% 'o','MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 2
figure(1)
plot([ones(1,length(w_percentage_pairs)).*(1+(rand(size(w_percentage_pairs))-0.4)/5)], (w_percentage_pairs),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on

plot([3*ones(1,length(w_percentage_pb))], (w_percentage_pb),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([2*ones(1,length(w_percentage_pre))], (w_percentage_pre),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([2*ones(1,length(w_percentage_pre));3* ones(1,length(w_percentage_pre))],[w_percentage_pre,w_percentage_pb]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on

plot([3*ones(1,length(w_percentage_pb))], (w_percentage_pb),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([4*ones(1,length(w_percentage_post))], (w_percentage_post),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([3*ones(1,length(w_percentage_pre));4* ones(1,length(w_percentage_pre))],[w_percentage_pb,w_percentage_post]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on

plot([0.75 1.25],[avg_w_percentage_pairs avg_w_percentage_pairs],'Color',[195 40 85]/256,'LineWidth',7)
plot([1.75 2.25],[avg_w_percentage_pre avg_w_percentage_pre],'Color',[195 40 85]/256,'LineWidth',7)
plot([2.75 3.25],[avg_w_percentage_pb avg_w_percentage_pb],'Color',[195 40 85]/256,'LineWidth',7)
plot([3.75 4.25],[avg_w_percentage_post avg_w_percentage_post],'Color',[195 40 85]/256,'LineWidth',7)

hold on
xlim([0 5])
ylim([0 40])
ylabel('% of whistle songs')
box off
xticks([1 2 3 4])
xticklabels({'Pairs','PRE','PB','POST'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on

%% Percentage of whistle songs: PRE-PB-POST
figure (2)
plot([2*ones(1,length(w_percentage_pb))], (w_percentage_pb),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([1*ones(1,length(w_percentage_pre))], (w_percentage_pre),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([1*ones(1,length(w_percentage_pre));2* ones(1,length(w_percentage_pre))],[w_percentage_pre,w_percentage_pb]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on

plot([2*ones(1,length(w_percentage_pb))], (w_percentage_pb),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([3*ones(1,length(w_percentage_post))], (w_percentage_post),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([2*ones(1,length(w_percentage_pre));3* ones(1,length(w_percentage_pre))],[w_percentage_pb,w_percentage_post]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on

plot([0.75 1.25],[avg_w_percentage_pre avg_w_percentage_pre],'Color',[195 40 85]/256,'LineWidth',7)
plot([1.75 2.25],[avg_w_percentage_pb avg_w_percentage_pb],'Color',[195 40 85]/256,'LineWidth',7)
plot([2.75 3.25],[avg_w_percentage_post avg_w_percentage_post],'Color',[195 40 85]/256,'LineWidth',7)

hold on
xlim([0 4])
ylim([0 50])
ylabel('% of whistle songs')
box off
xticks([1 2 3 ])
xticklabels({'PRE','PB','POST'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on
 set(gca,'TickDir','out')
%% percentage whistle matches over total songs (JUST WHISTLE SONGS!)

total_perc_whist_matches_pairs=(numbers_pairs(:,3)./numbers_pairs(:,1))*100;

%% percentage whistle matches over whistle songs (JUST WHISTLE SONGS!)

perc_whist_matches_whistle_pairs=(numbers_pairs(:,3)./numbers_pairs(:,2))*100;

%%
figure
plot([ones(1,length(total_perc_whist_matches_pairs)).*(1+(rand(size(total_perc_whist_matches_pairs))-0.4)/2)], (total_perc_whist_matches_pairs),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([0.75 1.25],[mean(total_perc_whist_matches_pairs) mean(total_perc_whist_matches_pairs)],'Color',[195 40 85]/256,'LineWidth',7)

xlim([0 2])
ylim([0 25])
ylabel('% of whistle matches')
xticks([1 ])
xticklabels({'Pairs'})
%legend(sprintf('p = %0.3f',p_anti))
%title('percentage of whistle matches over total songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on
box off
%%
figure
subplot(1,2,1)
plot([ones(1,length(w_percentage_pairs)).*(1+(rand(size(w_percentage_pairs))-0.4)/2)], (w_percentage_pairs),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[1 1 1], 'MarkerSize', 15)
hold on
plot([0.75 1.25],[avg_w_percentage_pairs avg_w_percentage_pairs],'Color',[195 40 85]/256,'LineWidth',7)
hold on
xlim([0 2])
ylim([0 35])
ylabel('% of whistle songs over total songs')
box off
xticks([1  ])
xticklabels({'Pairs'})
%legend(sprintf('p = %0.3f',p_anti))
%title('percentage of whistle songs')
set(gca,'linewidth',1,'FontSize', 14)
hold on
set(gca,'TickDir','out')


%figure(5)
subplot(1,2,2)
plot([ones(1,length(perc_whist_matches_whistle_pairs)).*(1+(rand(size(perc_whist_matches_whistle_pairs))-0.4)/1.5)], (perc_whist_matches_whistle_pairs),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[1 1 1], 'MarkerSize', 15)
hold on
plot([0.75 1.25],[mean(perc_whist_matches_whistle_pairs) mean(perc_whist_matches_whistle_pairs)],'Color',[195 40 85]/256,'LineWidth',7)

xlim([0 2])
ylim([0 100])
ylabel('% of whistle matches')
xticks([1 ])
xticklabels({'Pairs'})
%legend(sprintf('p = %0.3f',p_anti))
%title('percentage of whistle matches over whistle songs')
set(gca,'linewidth',1,'FontSize', 14)
hold on
set(gca,'TickDir','out')
box off
%%
[h,p]=ranksum(perc_whist_matches_whistle_pairs,w_match_perc)



figure
plot([ones(1,length(perc_whist_matches_whistle_pairs)).*(1+(rand(size(perc_whist_matches_whistle_pairs))-0.4)/5)], (perc_whist_matches_whistle_pairs),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on

plot([2*ones(1,length(w_match_perc)).*(1+(rand(size(w_match_perc))-0.6)/5)], (w_match_perc),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([0.75 1.25],[mean(perc_whist_matches_whistle_pairs) mean(perc_whist_matches_whistle_pairs)],'Color',[195 40 85]/256,'LineWidth',7)
plot([1.75 2.25],[mean(w_match_perc) mean(w_match_perc)],'Color',[195 40 85]/256,'LineWidth',7)

hold on
xlim([0 4])
ylim([0 100])
ylabel('whistle songs used for whistle song type matching (%)')
box off
xticks([1 2 3 ])
xticklabels({'Pairs','PB Germany', 'PB Gambia'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on
 set(gca,'TickDir','out')
%%
x=[31 33 2 14 7 3 2 10 2 3];
y=[10 15 2 7 3 0 2 0 0  0];

w=x+y;

per=(x./w)*100;

plot([2.75 3.25],[mean(per) mean(per)],'Color',[195 40 85]/256,'LineWidth',7)
plot([3*ones(1,length(per)).*(1+(rand(size(per))-0.4)/5)], (per),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)


[h,p]=ranksum(per,w_match_perc)


%%

[h,p]=ranksum(total_perc_whist_matches_pairs,w_match_perc_over_tot)


figure
plot([ones(1,length(total_perc_whist_matches_pairs)).*(1+(rand(size(total_perc_whist_matches_pairs))-0.4)/5)], (total_perc_whist_matches_pairs),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on

plot([2*ones(1,length(w_match_perc_over_tot)).*(1+(rand(size(w_match_perc_over_tot))-0.4)/5)], (w_match_perc_over_tot),...
    'o','MarkerFaceColor', [255 175 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([0.75 1.25],[mean(total_perc_whist_matches_pairs) mean(total_perc_whist_matches_pairs)],'Color',[195 40 85]/256,'LineWidth',7)
plot([1.75 2.25],[mean(w_match_perc_over_tot) mean(w_match_perc_over_tot)],'Color',[195 40 85]/256,'LineWidth',7)

hold on
xlim([0 3])
ylim([0 50])
ylabel('% of whistle matches over total songs')
box off
xticks([1 2 ])
xticklabels({'Pairs','PB'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on

%%
