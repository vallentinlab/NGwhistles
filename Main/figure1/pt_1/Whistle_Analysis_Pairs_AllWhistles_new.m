%% Whistle analysis 2.0 for PAIRS %%
clear all
close all
clc

% Read all bird's data...
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*ed.csv');
whistle_table = [];

bird_pitch = [];
bird_dur = [];
n_single_whistles=[];
percentile_birds=[];

delta_freq=[];
delta_dur=[];

collect_peaks_all= nan(25,25);
collect_delta=nan(25,25);

%% Load data/ adress folder %%
for bird=1:20
    
    
    
    
if bird == 1
    fprintf("Count with me! :)\n");
end
bird

filename= char(strcat(pathname,'\',filelist(bird,1).name));
[data,text] = xlsread(filename);

%% Read Excel %%

% % All whistles
bird_start = data(:,3);
bird_end = data(:,4);
bird_freq = data(:,5);
bird_duration = (bird_end - bird_start);
bird_dur = [bird_dur; bird_duration];
bird_pitch = [bird_pitch; bird_freq]; % accumulated pitch

% histogram how many single whistles

A=[data(:,2);1];
B=[];
for i=1:length(A)-1;
    if A(i+1)<=A(i); 
        B=[B; A(i)];
    end
end

n_single_whistles=[n_single_whistles;B];

P = prctile(bird_freq,[25 75],"all")';
percentile_birds=[percentile_birds; P];

% delta freq

C=[];
 for ii=1:height(data)-1;
if A(ii+1)>A(ii);
    C=[C; data(ii,5)-data(ii+1,5)];
end

 end
 
 delta_freq=[delta_freq;C];
 
 % delta duration
 D=[];
 for iii=1:height(data)-1;
if A(iii+1)>A(iii);
    D=[D; (data(iii,4)-(data(iii,3)))-(data(iii+1,4)-(data(iii+1,3)))];
end
 end
 
  delta_dur=[delta_dur;D];

%   figure (bird+2)
%   subplot(2,1,1)
%   histogram(bird_freq,300,...
%     'FaceColor', 'k', 'EdgeColor', [0 100 100]/256,'Normalization','probability')
%     xlim([0 10000]);
% hold on
% subplot(2,1,2)
% ksdensity(bird_freq, 'bandwidth', 50)
% [f, Xi, u] = ksdensity(bird_freq, 'bandwidth', 50);
%      xlim([0 10000]);
%     box off
%    [a,b]=findpeaks(f);
%         c=[a;b];
%         cc=sortrows(c');
%         freq1=Xi(cc(end,2));
%         hold on
%         for s=1:size(cc,1)
%         plot(Xi(cc(s,2)),cc(s,1), 'o','Color', 'r' )
%         end
%         collect_delta(bird,1:length(diff(Xi(cc(:,2)))))=diff(sort(Xi(cc(:,2))));
%         collect_peaks_all(bird,1:length((Xi(cc(2:end,2)))))=(sort(Xi(cc(2:end,2))));

  counts = hist(bird_freq,300);
   B =    smoothdata([ 0 0 0 0 0 0 0 0 0 0 counts/length(bird_freq) 0 0 0 0 0 0 0 0 0 0 ],'gaussian',15);
    xx=linspace(min(bird_freq),max(bird_freq),300);
    steps=xx(2)-xx(1);
    C=linspace(min(bird_freq)-10*steps, max(bird_freq)+10*steps,320);
    [a,b]=findpeaks(B);
        c=[a;b];
        cc=sortrows(c');
 
%  xlim([0 10000])
%      hold on
%         for s=1:size(cc,1)
%         plot(C(cc(s,2)),cc(s,1), 'o','Color', 'r' )
%         end
%  subplot(2,2,3);scatter((sort(C(cc(2:end,2)))),diff(sort(C(cc(:,2)))))
% xlim([0 10000])
% ylim([0 1500])
% subplot(2,2,4);plot((sort(C(cc(2:end,2)))),diff(sort(C(cc(:,2)))))
% ylim([0 1500])
% xlim([0 10000])
 collect_delta(bird,1:length(diff(sort(C(cc(:,2))))))=diff(sort(C(cc(:,2))));
 collect_peaks_all(bird,1:length(((C(cc(2:end,2))))))=(sort(C(cc(2:end,2))));
% %   
 

% %   
%   
%  
% figure (2)
% subplot(20,1,bird)
% ksdensity(bird_freq, 'bandwidth', 50)
% [f, Xi, u] = ksdensity(bird_freq, 'bandwidth', 50);
%      xlim([0 10000]);
%     box off
%    [a,b]=findpeaks(f);
%         c=[a;b];
%         cc=sortrows(c');
%         freq1=Xi(cc(end,2));
% 
%         hold on
%         for s=1:size(cc,1)
%         plot(Xi(cc(s,2)),cc(s,1), 'o','Color', 'r' )
%         end


 
end % for loop, all birds

%% Frequency distribution histogram

bird_pitch = bird_pitch(~isnan(bird_pitch));
length(bird_pitch)
[perc_20_80 ]=prctile(bird_pitch, [20 80]);
        
median(bird_pitch)/10

figure
subplot(2,1,2)
raincloud_plot_smooth(bird_pitch, 'line_width', 2,...
    'line_color', [35 70 20]/256, 'bxcl', [35 70 20]/256, 'color', [142 197 69]/256);
%     ylim([-0.0006 0.0006]);
%set(gca, 'XScale', 'log')
box off
 xlim([0 10000]);
     ylabel('probability');
xlabel('Pitch of whistle syllables (Hz)');
set(gca,'TickDir','out')
xline(perc_20_80(1), '--r')
xline(perc_20_80(2), '--r')

subplot(2,1,1)
histogram(bird_pitch,300,...
    'FaceColor', [0 100 100]/256, 'EdgeColor', [0 100 100]/256,'Normalization','probability')
    
    ylabel('probability');
      xlim([0 10000]);  
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
    box off
xline(perc_20_80(1), '--r')
xline(perc_20_80(2), '--r')
ylim([0 0.03]);
set(gca,'TickDir','out')

%%
% figure
% histogram(bird_pitch,300,'Normalization','probability','DisplayStyle','stairs')
% 
%  [f,xi] = ksdensity(bird_pitch); 
%  figure
% plot(xi,f);
%    figure
%    ecdf(bird_pitch)
%% Number of single whistles
 [h,p] = kstest((n_single_whistles))
 
 length(bird_pitch)
 n_single_whistles = n_single_whistles(~isnan(n_single_whistles));

 
 
figure
hold on
histogram(n_single_whistles,...
    'FaceColor', [0 100 100]/256, 'EdgeColor', [0 100 100]/256,'Normalization','probability')
    ylabel('No. of whistle songs');
    xlabel('No. of single whistles within whistle song');
    xlim([0 35])
    ylim([0 0.20])
    mean(n_single_whistles);
    box off
    hold on
xline(median(n_single_whistles), '--')
std(n_single_whistles)
set(gca,'TickDir','out')
%% Delta of frequency and duration

% Frequency
figure
subplot(1,2,1)
% violinplot(delta_freq)
[p,h]= ranksum(0,delta_freq)
[p,h]= ttest(delta_freq)

%     ylabel('Delta of pitch between single whistles (Hz)');
%     %ylim([-1000 1000])
%     xlim([0.5 1.5])
%     axis square
%     set(gca,'xtick',[]);
%     set(gca,'xticklabel',[]);
    
raincloud_plot_smooth(delta_freq, 'line_width', 2,...
    'line_color', [35 70 20]/256, 'box_on', 0, 'bxcl', [35 70 20]/256, 'color', [142 197 69]/256);
    ylim([-0.08 0.08]);
     xlim([-1500 1500]);
     xlabel('Delta of frequency between single whistles (Hz)');
     % set(gca,'yticklabel',[{''};{''};{'0'};{'0.005'};{'0.01'}]);
    axis square
    box off
xline(0, '--')
set(gca,'TickDir','out')
    
std(delta_freq, 'omitnan')
mean(delta_freq, 'omitnan')
%% Duration
hold on
subplot(1,2,2)
% violinplot(delta_dur)
% ranksum(0,delta_dur)
%     ylabel('Delta of duration between single whistles (s)');
%     ylim([-1 1])
%     xlim([0.5 1.5])
%     axis square
%     set(gca,'xtick',[]);
%     set(gca,'xticklabel',[]);

raincloud_plot_smooth(delta_dur, 'line_width', 2,...
    'line_color', [35 70 20]/256,'box_on', 0, 'bxcl', [35 70 20]/256, 'color', [142 197 69]/256);
     ylim([-0.08 0.08]);
     xlim([-0.8 0.8]);
     xlabel('Delta of duration between single whistles (s)');
     %set(gca,'yticklabel',[{''};{''};{''};{''};{'0'};{'5'};{'10'};{'15'};{'20'}]);
    axis square
    box off
xline(0, '--')
set(gca,'TickDir','out')

std(delta_dur, 'omitnan')
mean(delta_dur, 'omitnan')
% 
% find(delta_freq>1000)
% 
% filename
% [x,z]=find(C<-1000)
% [x,z]=find(C>1000)
% [x,y]=find(bird_pitch>9000)
%%


REcollect_peaks_all=reshape(collect_peaks_all',[625,1]);
REcollect_delta=(reshape(collect_delta',[625,1]));

REcollect_peaks_all=REcollect_peaks_all(~isnan(REcollect_peaks_all));
        REcollect_delta=REcollect_delta(~isnan(REcollect_delta));  

[rho,pval] = corr(REcollect_delta,REcollect_peaks_all,'Type','Pearson','Rows','complete')

figure;scatter(REcollect_peaks_all,REcollect_delta)
 yline(mean(REcollect_delta,'omitnan'), '--')
   yline( median(REcollect_delta,'omitnan'), ':')
set(gca, 'YScale', 'log')

    tpb = table(REcollect_peaks_all,...
        REcollect_delta);
mpb = fitlm(tpb,'REcollect_delta ~ REcollect_peaks_all')

figure
ppb = plot(mpb,'color','k',...
    'Marker','o','MarkerFaceColor', [102 162 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(ppb(2),'Color',[195 40 85]/256,'LineWidth',5);
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    axis square
    
    
    figure
    histogram(REcollect_delta,50,'Normalization','probability')
    box off
    axis square
set(gca, 'XScale', 'log')

    mean(REcollect_delta,'omitnan')
    median(REcollect_delta,'omitnan')
    
    
    DagosPtest(log(REcollect_delta))
    
    
    %% preferred range
    REcollect_peaks_all=REcollect_peaks_all(~isnan(REcollect_peaks_all));
        REcollect_delta=REcollect_delta(~isnan(REcollect_delta));  
  RE=[REcollect_peaks_all REcollect_delta];
    
pref_delta=[];
pref_peak=[];
for ww=1:size(RE,1)
    if RE(ww,1) < perc_20_80(2) && RE(ww,1)>perc_20_80(1)
        pref_delta=[pref_delta; RE(ww,2)];
        pref_peak=[pref_peak; RE(ww,1)];
    end
end

 tpb = table(pref_peak,...
        pref_delta);
mpb = fitlm(tpb,'pref_delta ~ pref_peak')

figure
ppb = plot(mpb,'color','k',...
    'Marker','o','MarkerFaceColor', [102 162 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(ppb(2),'Color',[195 40 85]/256,'LineWidth',5);
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    axis square
    %%
    
    figure;histogram(log(bird_pitch),300)
    counts = hist((bird_pitch),300)
     
   B =    smoothdata([ 0 0 0 0 0 0 0 0 0 0 counts/length(bird_pitch) 0 0 0 0 0 0 0 0 0 0 ],'gaussian',15);
    xx=linspace(min(bird_pitch),max(bird_pitch),300);
    steps=xx(2)-xx(1);

    C=linspace(min(bird_pitch)-10*steps, max(bird_pitch)+10*steps,320);
    
    
    
    [acf,lags] = autocorr(counts,299);
    figure;plot(acf);
    [aa,bb]=findpeaks(acf);
    hold on;
    plot(bb,aa, 'o','Color', 'm' );
    
    
    [a,b]=findpeaks(B);
        c=[a;b];
        cc=sortrows(c')
    
    figure
    subplot(2,1,2)
  
raincloud_plot_smooth(bird_pitch, 'line_width', 2,...
    'line_color', [35 70 20]/256, 'bxcl', [35 70 20]/256, 'color', [142 197 69]/256);
ylim([-0.02 0.02]);
set(gca,'TickDir','out')
box off
   for s=1:size(cc,1)
        plot(C(cc(s,2)),cc(s,1)+0.001, '*','Color', 'r' )
        end
 xlim([0 10000]);
    % xline(perc_20_80(1), '--r')
%xline(perc_20_80(2), '--r')
 
subplot(2,1,1)
histogram(bird_pitch,300,...
    'FaceColor', [0 100 100]/256, 'EdgeColor', [0 100 100]/256,'Normalization','probability')
    xlabel('Whistle frequency (Hz)');
    ylabel('No. of single whistles');
      xlim([0 10000]);  
ylim([0 0.03]);
set(gca,'TickDir','out')
    box off
 hold on
   % xline(perc_20_80(1), '--r')
%xline(perc_20_80(2), '--r')
    
   %%
%     subplot(3,1,3)
%  plot(C, B)
%  xlim([1000 10000])
%      hold on
%         for s=1:size(cc,1)
%         plot(C(cc(s,2)),cc(s,1), 'o','Color', 'r' )
%         end
%     xline(perc_20_80(1), '--r')
% xline(perc_20_80(2), '--r')
% box off
%     yline(0.005)    
        
        
DDelta =(diff(sort(C(cc(:,2)))))'
 PPeak=(sort(C(cc(2:end,2))))'
 tpb = table(PPeak,...
        DDelta);
mpb = fitlm(tpb,'DDelta~ PPeak ')

figure
ppb = plot(mpb,'color','k',...
    'Marker','o','MarkerFaceColor', [102 162 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(ppb(2),'Color',[195 40 85]/256,'LineWidth',5);
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    axis square
    
    
figure;scatter(PPeak,DDelta)
yline(mean(DDelta),'--')

figure
scatterhist(PPeak,DDelta,'Kernel','on',...
    'Direction','out','Marker','*','MarkerSize',10);
yline(mean(DDelta),'--')
xlim([0 10000])
ylim([0 1000])
box off
axis square
set(gca,'TickDir','out')


figure;histogram(DDelta,7)


    DagosPtest(log(DDelta))
%%
sil=[];
perc_sil=[];

 for r=3:30
    
     [idx,C] = kmeans(bird_pitch,r,'Replicates',100);

     
     [s,h]=silhouette(bird_pitch,idx);
     
sil=[sil,     mean(s)];
     

[w,ww]=find(s<0);

perc_sil=[ perc_sil sum(ww)/length(s)*100];

     
 end

sil=[nan nan sil]

figure
plot(sil)
hold on
plot(10,max(sil), 'or')
xlim([1 35])
ylim([0.7 0.8])
box off
axis square
set(gca,'TickDir','out')
  xlabel('Number of Clusters');
    ylabel('Silhouette criterion');
    
%%

% D = pdist(bird_pitch);
% 
% Z = linkage(D,'centroid');
% dendrogram(Z,0)
% cutoff = 500
% dendrogram(Z,'ColorThreshold',cutoff)
% 
%%
% sil=[];
% perc_sil=[];
% 
%  for r=3:30
%     
% T = clusterdata(bird_pitch,'Linkage','average','Cutoff');
% 
%      
%      [s,h]=silhouette(bird_pitch,T);
%      
% sil=[sil,     mean(s)];
%      
% 
% [w,ww]=find(s<0);
% 
% perc_sil=[ perc_sil sum(ww)/length(s)*100];
% 
%      
%  end
% 
% sil


%%


[idx,C] = kmeans(bird_pitch,10,'Replicates',100);

figure
[s,h]=silhouette(bird_pitch,idx)

  counts = hist(bird_pitch,300);
   B =    smoothdata([ 0 0 0 0 0 0 0 0 0 0 counts/length(bird_pitch) 0 0 0 0 0 0 0 0 0 0 ],'gaussian',15);
    xx=linspace(min(bird_pitch),max(bird_pitch),300);
    steps=xx(2)-xx(1);
    D=linspace(min(bird_pitch)-10*steps, max(bird_pitch)+10*steps,320);



figure
subplot(2,1,1)
plot(D, B  , 'LineWidth', 3); 
hold on
histogram(bird_pitch,300,...
    'FaceColor', [0 100 100]/256, 'EdgeColor', [0 100 100]/256,'Normalization','probability')
    ylabel('probability');     
set(gca, 'XScale', 'log')
    box off
xline(perc_20_80(1), '--r')
xline(perc_20_80(2), '--r')
ylim([0 0.03]);
xlim([1000 10000]); 
set(gca,'TickDir','out')
hold on
for ii=1:length(C)
    plot(C(ii),0.025, '*')
hold on
end

subplot(2,1,2)
for ii=1:length(C)
    oo=[rand rand rand];
    scatter(bird_pitch(idx==ii),(10+(rand(length((1:length(bird_pitch(idx==ii)))),1))*30),...
        'MarkerEdgeColor', 'none','MarkerFaceColor',oo);
hold on
end
hold on
for ii=1:length(C)
    plot(C(ii),45, '*')
hold on
end
xline(perc_20_80(1), '--r')
xline(perc_20_80(2), '--r')

set(gca,'TickDir','out')
set(gca, 'xScale', 'log')
 ylim([0  50 ]); 
 xlim([1000 10000]); 

 
 [h,p]=kstest(bird_pitch)
 
 %% 
 
 low_limit=[];
 
 high_limit=[];
 
 for ii=1:length(C)
     
     low_limit=[ low_limit;    min(bird_pitch(idx==ii))];
    high_limit=[high_limit; max(bird_pitch(idx==ii))];

end
hold on
 
limits=[low_limit, high_limit]
 
 
figure
 scatter(bird_pitch,(10+(rand(length((1:length(bird_pitch))),1))*1),...
        'MarkerEdgeColor', 'none','MarkerFaceColor',oo);

    for t=1:length(low_limit)
        xline(low_limit(t))
    hold on
    end
    
 %% 
 
 D=sort(C);
 [h,p] = lillietest(log(diff(D)))

 

 DagosPtest(log(diff(D)))

 [f, Xi, u] = ksdensity(diff(D),'Support','positive')
 
 
figure
subplot(1,2,1)
e=[rand rand rand]
scatter(D(2:end),(diff(sort(C))),'Marker','*',...
     'MarkerFaceColor',[e])
ylim([100 1000])
xlim([1000 10000])
box off
axis square
set(gca,'TickDir','out')
 ylim([100 2500]);
%set(gca, 'yScale', 'log')
 %set(gca, 'xScale', 'log')
  xlabel('Pitch of cluster i (Hz)');
    ylabel('Delta Pitch (cluster i - cluster i-1) (Hz) ');
    
 
 
figure
scatter(D(1:end-1),D(2:end),200,'Marker','*',...
     'MarkerFaceColor',[e])
 ylim([1000 10000])
xlim([1000 10000])
box off
axis square
set(gca,'TickDir','out')
%  set(gca, 'yScale', 'log')
%   set(gca, 'xScale', 'log')
  xlabel('Pitch of cluster i (Hz)');
    ylabel('Delta Pitch (cluster ');
    
    
    %%


figure
subplot(1,2,1)
scatter(D(2:end),diff(sort(C)))
ylim([0 1000])
xlim([0 10000])
box off
axis square
set(gca,'TickDir','out')
subplot(1,2,2)
violinplot((diff(sort(C))))
ylim([0 1000])
xlim([0 2])
box off
axis square
set(gca,'TickDir','out')