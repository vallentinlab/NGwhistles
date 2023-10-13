%% Whistle analysis 2.0 %%
clear all
close all
clc

% Set 1:
pb1_stim = []; pb1_bird = []; pb1_dur = [];
post1_stim = []; post1_bird = []; post1_dur = [];
pre1_stim = []; pre1_bird = []; pre1_dur = []; tempo1=[];wgap1=[];
% Set 2:
pb2_stim = []; pb2_bird = []; pb2_dur = [];
post2_stim = []; post2_bird = []; post2_dur = [];
pre2_stim = []; pre2_bird = []; pre2_dur = []; tempo2=[];wgap2=[];

% Read all bird's data...
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*.xlsx');
whistle_table = [];

%% Load data/ adress folder %%
for bird=1:11

if bird == 1
    fprintf("Count with me! :)\n");
end
bird

filename= char(strcat(pathname,'\',filelist(bird,1).name));
[~,sheet_name]=xlsfinfo(filename);
for k=1:numel(sheet_name)
  data{k}=xlsread(filename,sheet_name{k});
%   if isempty(data{1,k})==1
%       data{1,k} = nan;
%   end
end
% Example: x = data{1,2}(:,3);
% data{1,2} will load the 2nd sheet of the 1st file in filelist
% ...(:,3) will get all the rows from the 3rd column of this sheet

%% Read Excel %%
% Damn, this is long. There should be an easier way to do this
% Set 1
% Playback: pb10, pb11, pb12, pb13
if isempty(data{1,1})==1
    pb10_pitch = nan;    pb10_freq = nan;
    pb10_start = nan;    pb10_end = nan;
else
    pb10_pitch = data{1,1}(:,3);    pb10_freq = data{1,1}(:,14);
    pb10_start = data{1,1}(:,1);    pb10_end = data{1,1}(:,2);
end
if isempty(data{1,2})==1
    pb11_pitch = nan;    pb11_freq = nan;
    pb11_start = nan;    pb11_end = nan;
else
    pb11_pitch = data{1,2}(:,3);    pb11_freq = data{1,2}(:,14);
    pb11_start = data{1,2}(:,1);    pb11_end = data{1,2}(:,2);
end
if isempty(data{1,3})==1
    pb12_pitch = nan;    pb12_freq = nan;
    pb12_start = nan;    pb12_end = nan;
else
    pb12_pitch = data{1,3}(:,3);    pb12_freq = data{1,3}(:,14);
    pb12_start = data{1,3}(:,1);    pb12_end = data{1,3}(:,2);
end
if isempty(data{1,4})==1
    pb13_pitch = nan;    pb13_freq = nan;
    pb13_start = nan;    pb13_end = nan;
else
    pb13_pitch = data{1,4}(:,3);    pb13_freq = data{1,4}(:,14);
    pb13_start = data{1,4}(:,1);    pb13_end = data{1,4}(:,2);
end

% Postcontrol: post1
if isempty(data{1,5})==1
    post1_pitch = nan;    post1_freq = nan;
    post1_start = nan;    post1_end = nan;
else
    post1_pitch = data{1,5}(:,3);    post1_freq = data{1,5}(:,14);
    post1_start = data{1,5}(:,1);    post1_end = data{1,5}(:,2);
end

% Precontrol: pre1
if isempty(data{1,6})==1
    pre1_pitch = nan;    pre1_freq = nan;
    pre1_start = nan;    pre1_end = nan;
else
    pre1_pitch = data{1,6}(:,3);    pre1_freq = data{1,6}(:,14);
    pre1_start = data{1,6}(:,1);    pre1_end = data{1,6}(:,2);
end
 
% Set 2
% Playback: pb20, pb21, pb22
if isempty(data{1,7})==1
    pb20_pitch = nan;    pb20_freq = nan;
    pb20_start = nan;    pb20_end = nan;
else
    pb20_pitch = data{1,7}(:,3);    pb20_freq = data{1,7}(:,14);
    pb20_start = data{1,7}(:,1);    pb20_end = data{1,7}(:,2);
end
if isempty(data{1,8})==1
    pb21_pitch = nan;    pb21_freq = nan;
    pb21_start = nan;    pb21_end = nan;
else
    pb21_pitch = data{1,8}(:,3);    pb21_freq = data{1,8}(:,14);
    pb21_start = data{1,8}(:,1);    pb21_end = data{1,8}(:,2);
end
if isempty(data{1,9})==1
    pb22_pitch = nan;    pb22_freq = nan;
    pb22_start = nan;    pb22_end = nan;
else
    pb22_pitch = data{1,9}(:,3);    pb22_freq = data{1,9}(:,14);
    pb22_start = data{1,9}(:,1);    pb22_end = data{1,9}(:,2);
end

% Postcontrol: post2
if isempty(data{1,10})==1
    post2_pitch = nan;    post2_freq = nan;
    post2_start = nan;    post2_end = nan;
else
    post2_pitch = data{1,10}(:,3);    post2_freq = data{1,10}(:,14);
    post2_start = data{1,10}(:,1);    post2_end = data{1,10}(:,2);
end

% Precontrol: pre20, pre21
if isempty(data{1,11})==1
    pre20_pitch = nan;    pre20_freq = nan;
    pre20_start = nan;    pre20_end = nan;
else
    pre20_pitch = data{1,11}(:,3);    pre20_freq = data{1,11}(:,14);
    pre20_start = data{1,11}(:,1);    pre20_end = data{1,11}(:,2);
end
if isempty(data{1,12})==1
    pre21_pitch = nan;    pre21_freq = nan;
    pre21_start = nan;    pre21_end = nan;
else
    pre21_pitch = data{1,12}(:,3);    pre21_freq = data{1,12}(:,14);
    pre21_start = data{1,12}(:,1);    pre21_end = data{1,12}(:,2);
end

% Storing data from all birds in new variables
pb1_stim = [pb1_stim; pb10_freq; pb11_freq; pb12_freq; pb13_freq];
pb1_bird = [pb1_bird; pb10_pitch; pb11_pitch; pb12_pitch; pb13_pitch];
pb1_dur = [pb1_dur; (pb10_end - pb10_start); (pb11_end - pb11_start); (pb12_end - pb12_start); (pb13_end - pb13_start)];
pb2_stim = [pb2_stim; pb20_freq; pb21_freq; pb22_freq];
pb2_bird = [pb2_bird; pb20_pitch; pb21_pitch; pb22_pitch];
pb2_dur = [pb2_dur; (pb20_end - pb20_start); (pb21_end - pb21_start); (pb22_end - pb22_start)];
post1_stim = [post1_stim; post1_freq];
post1_bird = [post1_bird; post1_pitch];
post1_dur = [post1_dur; (post1_end - post1_start)];
post2_stim = [post2_stim; post2_freq];
post2_bird = [post2_bird; post2_pitch];
post2_dur = [post2_dur; (post2_end - post2_start)];
pre1_stim = [pre1_stim; pre1_freq];
pre1_bird = [pre1_bird; pre1_pitch];
pre1_dur = [pre1_dur; (pre1_end - pre1_start)];
tempo1=[tempo1; (pre1_start(2:end) - pre1_start(1:end-1))];
wgap1=[wgap1; (pre1_end(1:end-1) - pre1_start(2:end))];

pre2_stim = [pre2_stim; pre20_freq; pre21_freq];
pre2_bird = [pre2_bird; pre20_pitch; pre21_pitch];
pre2_dur = [pre2_dur; (pre20_end - pre20_start); (pre21_end - pre21_start)];
tempo2=[tempo2; (pre20_start(2:end)- pre20_start(1:end-1));(pre21_start(2:end) - pre21_start(1:end-1))];
wgap2=[wgap2; (pre20_end(1:end-1) - pre20_start(2:end));(pre21_end(1:end-1) - pre21_start(2:end))];

end % for loop, all birds

% Save table
% [~,name,~]=fileparts(pwd); % obtains name of the folder
% Warning: If folder name contains dots, Matlab will not that read from
% that point on
% tablename = [name,'_Pitch_Table.csv'];
% writetable(whistle_table, tablename, 'WriteRowNames', true);

%% Visualizing data %%

post_bird=[post1_bird; post2_bird];
pre_bird=[pre1_bird; pre2_bird];
pb_stim=[pb1_stim; pb2_stim];
pb_bird=[pb1_bird;pb2_bird];
%%

length(post_bird)

%% Stimuli VS bird frequency (all) %%

pb_bird=[pb1_bird; pb2_bird];
pb_stim=[pb1_stim; pb2_stim];

tpb = table(pb_bird, pb_stim);
mpb = fitlm(tpb,'pb_bird ~ pb_stim')

figure(5)
ppb = plot(mpb,'color','r','Marker','o','MarkerSize',4,'MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','r','LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
    xlim([0 9000]);
    ylim([0 9000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')

[rho,pval] = corr(pb_bird, pb_stim,'Type','Pearson','Rows','complete')



%%


pb_bird=[pb1_bird; pb2_bird];
pb_stim=[pb1_stim; pb2_stim];

difff=abs(pb_stim-pb_bird)
tpb = table(difff, pb_stim);
mpb = fitlm(tpb,'difff ~ pb_stim')

figure
ppb = plot(mpb,'color','r','Marker','o','MarkerSize',4,'MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','r','LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
   
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')

[rho,pval] = corr(pb_bird, pb_stim,'Type','Pearson','Rows','complete')

%%
J=customcolormap([0 0.5 1], [1 1  0 ; 0 0 1;1 1 1 ]);
J=customcolormap([0 0.5 1], [0 0.2 0; 0 0.5 0;1 1 1 ]);

figure
 [values, centers] = hist3([   pb_stim difff],[100 100]);
imagesc(centers{:}, values.')

colorbar
axis xy 
% xlim([0 8000]);
%     ylim([0 8000]);
        axis square
colormap(J)
hold on
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
     legend off
%     xlim([0 9000]);
%     ylim([0 9000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')
 
 
xline(1735,'-') 
 xline(3652,'-')

%   yline(1735,'--')
%  yline(3652,'--')

% q=[0,9000];
% qq=[0,9000];
% hold on
%   plot (q,qq, ':k')


%% (Histogram) PRE-PB-POST freq. distribution %%
figure(1)
subplot(3,2,1)
violinplot(pre_bird);
    %title('PRE');
    set(gca,'xtick',[],'xticklabel',[]);
    xlim([0.5 1.5]);
    axis square

% [f, Xi, u] = ksdensity(pre_bird);
%           %[a,b]=max(f)
%           [a,b]=findpeaks(f);
%          c=[a;b];
%          cc=sortrows(c');
%      freq1=Xi(cc(end,2));
%      freq2=Xi(cc(end-1,2));
%       freq3=Xi(cc(end-2,2));
%      xline(freq1,'-','LineWidth',8,'Color', 'r' )
%                xline(freq2,'-','LineWidth',8,'Color', 'm' )
%                 xline(freq3,'-','LineWidth',8,'Color', 'k' )

subplot(3,2,2)
histogram(pre_bird,300,...
    'FaceColor', [80 160 50]/256, 'EdgeColor', [80 160 50]/256,'Normalization','probability');
    %title('PRE');
    xlim([0 9000]);
    ylim([0 0.1]);
    ylabel('No. of whistles');
    xlabel('Frequency (Hz)');

subplot(3,2,3)
histogram(pb_stim, 70);
    title('PB stimuli');

subplot(3,2,4)
histogram(pb_bird,300,...
    'FaceColor', [80 160 50]/256, 'EdgeColor', [80 160 50]/256,'Normalization','probability');
    title('Bird responses to PB');
    xlim([0 9000]);
    ylim([0 0.1]);
    ylabel('No. of whistles');
    xlabel('Frequency (Hz)');

subplot(3,2,5)
violinplot(post_bird);
    %title('POST');
    set(gca,'xtick',[],'xticklabel',[]);
    xlim([0.5 1.5]);
    axis square

subplot(3,2,6)
histogram(post_bird,300,...
    'FaceColor', [80 160 50]/256, 'EdgeColor', [80 160 50]/256,'Normalization','probability');
    %title('POST');
    xlim([0 9000]);
    ylim([0 0.1]);
    ylabel('No. of whistles');
    xlabel('Frequency (Hz)');

[p,h]=ranksum(pb_bird,pre_bird);
[p,h]=ranksum(pb_bird,post_bird);
[p,h]=ranksum(post_bird,pre_bird);

nanmedian(post_bird)
% figure
% violinplot(pre_bird)

%% (Raincloud) PRE-PB-POST freq. distribution %%

figure
subplot(3,1,1)
raincloud_plot(pre_bird, 'line_width', 2,...
    'cloud_edge_col', [35 70 20]/256, 'bxcl', [35 70 20]/256, 'color', [142 197 69]/256);
    %title('PRE');
    ylim([-0.0006 0.0006]);
    box off

     pb_bird=pb_bird(~isnan(pb_bird));
    
subplot(3,1,2)
raincloud_plot(pb_bird, 'line_width', 2,...
    'cloud_edge_col', [35 70 20]/256, 'bxcl', [35 70 20]/256, 'color', [142 197 69]/256);
    %title('Bird responses to PB');
    ylim([-0.0006 0.0006]);
    ylabel('No. of whistles');
    box off

subplot(3,1,3)
raincloud_plot(post_bird, 'line_width', 2,...
    'cloud_edge_col', [35 70 20]/256, 'bxcl', [35 70 20]/256, 'color', [142 197 69]/256);
    %title('POST');
    ylim([-0.0006 0.0006]);
    xlabel('Frequency (Hz)');
    box off
%% (Raincloud) Separated version

% PRE
figure(3)
subplot(2,1,1)
histogram(pre_bird,300,'Normalization','probability');
    title('PRE');
    xlim([0 9000]);
    %ylim([0 150]);
    ylabel('No. of whistles');
    box off
      
subplot(2,1,2)
raincloud_plot(pre_bird);
    title('PRE');
    ylim([-0.0006 0.0006]);
    xlim([0 9000]);
    xlabel('Frequency (Hz)');
    box off

% POST
figure(4)
subplot(2,1,1)
histogram(post_bird,400);
    title('POST');
    xlim([0 9000]);
    ylim([0 150]);
    ylabel('No. of whistles');
    box off
        
subplot(2,1,2)
raincloud_plot(post_bird);
    title('POST');
    ylim([-0.0006 0.0006]);
    xlim([0 9000]);
    xlabel('Frequency (Hz)');
    box off
   
%% Stimuli VS bird frequency (all) %%

pb_bird=[pb1_bird; pb2_bird];
pb_stim=[pb1_stim; pb2_stim];

tpb = table(pb_bird, pb_stim);
mpb = fitlm(tpb,'pb_bird ~ pb_stim')

figure(5)
ppb = plot(mpb,'color','r','Marker','o','MarkerSize',4,'MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','r','LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
    xlim([0 9000]);
    ylim([0 9000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')

[rho,pval] = corr(pb_bird, pb_stim,'Type','Pearson','Rows','complete')

%%
J=customcolormap([0 0.5 1], [1 1  0 ; 0 0 1;1 1 1 ]);
J=customcolormap([0 0.5 1], [0 0.2 0; 0 0.5 0;1 1 1 ]);

figure
 [values, centers] = hist3([  pb_stim pb_bird],[100 100]);
imagesc(centers{:}, values.')

colorbar

axis xy 
xlim([0 8000]);
    ylim([0 8000]);
        axis square
colormap(J)
hold on
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
    xlim([0 9000]);
    ylim([0 9000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')
  xline(1700,'--')
 xline(3700,'--')

q=[0,9000];
qq=[0,9000];
hold on
  plot (q,qq, ':k')
    %% divide in preferred and not
    
bird_pb=[pb_bird, pb_stim];


pref_bird=[];
pref_stim=[];
for ww=1:size(bird_pb,1)
    if bird_pb(ww,1) < 3652 && bird_pb(ww,1)>1735
        pref_bird=[pref_bird; bird_pb(ww,1)];
        pref_stim=[pref_stim; bird_pb(ww,2)];
    end
end

no_pref_bird=[];
no_pref_stim=[];
for ww=1:size(bird_pb,1)
    if bird_pb(ww,1) > 3652
        no_pref_bird=[no_pref_bird; bird_pb(ww,1)];
        no_pref_stim=[no_pref_stim; bird_pb(ww,2)];
        
    elseif bird_pb(ww,1)<1735
        no_pref_bird=[no_pref_bird; bird_pb(ww,1)];
        no_pref_stim=[no_pref_stim; bird_pb(ww,2)];
    end
end


no_pref_bird_pb=[no_pref_bird, no_pref_stim];
no_PREF_pref_bird=[];
no_PREF_pref_stim=[];
for ww=1:size(no_pref_bird_pb,1)
    if no_pref_bird_pb(ww,2) > 3652
        no_PREF_pref_bird=[no_PREF_pref_bird; no_pref_bird_pb(ww,1)];
        no_PREF_pref_stim=[no_PREF_pref_stim; no_pref_bird_pb(ww,2)];
        
    elseif no_pref_bird_pb(ww,2) <1735
        no_PREF_pref_bird=[no_PREF_pref_bird; no_pref_bird_pb(ww,1)];
        no_PREF_pref_stim=[no_PREF_pref_stim; no_pref_bird_pb(ww,2)];
    end
end


pref_bird_pb=[pref_bird, pref_stim]
PREF_pref_bird=[];
PREF_pref_stim=[];
for ww=1:size(pref_bird_pb,1)
    if pref_bird_pb(ww,2) > 3652
        PREF_pref_bird=[PREF_pref_bird; pref_bird_pb(ww,1)];
        PREF_pref_stim=[PREF_pref_stim; pref_bird_pb(ww,2)];
        
    elseif pref_bird_pb(ww,2) <1735
        PREF_pref_bird=[PREF_pref_bird; pref_bird_pb(ww,1)];
        PREF_pref_stim=[PREF_pref_stim; pref_bird_pb(ww,2)];
    end
end



%%
tpb = table(no_PREF_pref_bird, no_PREF_pref_stim);
mpb = fitlm(tpb,'no_PREF_pref_bird ~ no_PREF_pref_stim')



% tpb = table(PREF_pref_bird, PREF_pref_stim);
% mpb = fitlm(tpb,'PREF_pref_bird ~ PREF_pref_stim')

figure
hold on
ppb = plot(mpb,'color','r','Marker','o','MarkerFaceColor', 'b', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','b','LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
    xlim([0 9000]);
    ylim([0 9000]);
    title ('Freq. bird VS stimuli');
     axis square
    hold on
    
  

    %%
    
    idd=find(pb_stim>1700);
iddd=find(pb_stim<3700)
idddd=intersect(idd,iddd)

new_stim_pitch=pb_stim(idddd)
new_bird_pitch=pb_bird(idddd)

tpb = table(new_bird_pitch, new_stim_pitch);
mpbb = fitlm(tpb,'new_bird_pitch ~ new_stim_pitch')
figure(6)
ppbb = plot(mpbb,'color','k','Marker','o','MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppbb(2), 'Color','r','LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
    xlim([0 9000]);
    ylim([0 9000]);
    title ('Freq. bird VS stimuli');

%% Density plot (all) - violin version %%
post_bird = [post1_bird; post2_bird];
pre_bird = [pre1_bird; pre2_bird];

x = ones(length(pre_bird),1).*(1+(rand(length(pre_bird),1)-0.4)/5);
y = ones(length(post_bird),1).*(1+(rand(length(post_bird),1)-0.4)/5);
jj = [600:200:8600];

data_matrix = nan(10000);
data_matrix(1:length(pre_bird),1) = pre_bird;
data_matrix(1:length(post_bird),41) = post_bird;
 for ww=2:length (jj)
  ids=find(pb_stim==jj(ww));
  data_matrix(1:length(ids),ww)=pb_bird(ids);
  end
 
% figure(6)
% violinplot(data_matrix);
%     ylim([0 9000]);
%     hold on
%     xticks([1:41]);
% xticklabels({'PRE',convertCharsToStrings(jj(end-39)),convertCharsToStrings(jj(end-38)),...
%     convertCharsToStrings(jj(end-37)),convertCharsToStrings(jj(end-36)),convertCharsToStrings(jj(end-35)),...
%     convertCharsToStrings(jj(end-34)),convertCharsToStrings(jj(end-33)),convertCharsToStrings(jj(end-32)),...
%     convertCharsToStrings(jj(end-31)),convertCharsToStrings(jj(end-30)),convertCharsToStrings(jj(end-29)),...
%     convertCharsToStrings(jj(end-28)),convertCharsToStrings(jj(end-27)),convertCharsToStrings(jj(end-26)),...
%     convertCharsToStrings(jj(end-25)),convertCharsToStrings(jj(end-24)),convertCharsToStrings(jj(end-23)),...
%     convertCharsToStrings(jj(end-22)),convertCharsToStrings(jj(end-21)),convertCharsToStrings(jj(end-20)),...
%      convertCharsToStrings(jj(end-19)),convertCharsToStrings(jj(end-18)),convertCharsToStrings(jj(end-17)),...
%     convertCharsToStrings(jj(end-16)),convertCharsToStrings(jj(end-15)),convertCharsToStrings(jj(end-14)),...
%     convertCharsToStrings(jj(end-13)),convertCharsToStrings(jj(end-12)),convertCharsToStrings(jj(end-11)),...
%     convertCharsToStrings(jj(end-10)),convertCharsToStrings(jj(end-9)),convertCharsToStrings(jj(end-8)),...
%     convertCharsToStrings(jj(end-7)),convertCharsToStrings(jj(end-6)),convertCharsToStrings(jj(end-5)),...
%     convertCharsToStrings(jj(end-4)),convertCharsToStrings(jj(end-3)),convertCharsToStrings(jj(end-2)),...
%     convertCharsToStrings(jj(end)),'POST'}');
%     ylabel('Bird freq (Hz)');
%     xlabel('Stimulus freq (Hz)');

%% Density plot (all) - scatter version %%

figure(7)
% Precontrol
for qq=1
    subplot(1,length(jj),qq);
    s1 = scatter_kde(x+(10+(rand(length(pre_bird),1))*30),pre_bird);
    s1.SizeData = 15;
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('PRE');
    hold on
    [gg,bin]=hist(pre_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' );
    [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',7,'Color', 'r' ); % first peak
end

peak_all_bird=[];
peak_all_pb=[];
% Playback stimuli
for qq=2:length(jj);
    subplot(1,length(jj),qq);
    ids=find(pb_stim==jj(qq));
    pb_bird(ids);
    if (sum(pb_stim(ids))>0)==1;
        s2 = scatter_kde(pb_stim(ids)+(10+(rand(length(ids),1))*30),pb_bird(ids));
        s2.SizeData = 15;
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(jj(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',7); % pb frequency
      % yline(median(pb_bird(ids),'omitnan'),'-','Color', [195 40 85]/256, 'LineWidth',7);
      % yline(mean(pb_bird(ids),'omitnan'),'-','Color', [120 20 230]/230, 'LineWidth',7);

        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
      % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
        hold on
        [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
      % [a,b]=max(f)
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        %freq2=Xi(cc(end-1,2));
        yline(freq1,'-','LineWidth',7,'Color', 'r' ) % first peak
        % yline(freq2,'-','LineWidth',7,'Color', 'm' ) % second peak
peak_all_pb=[peak_all_pb;jj(qq)];
       % ranksum(jj(qq),pb_bird(ids));
        peak_all_bird=[peak_all_bird;freq1];
    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(jj(qq)));
        hold on
    end
hold on
end

% Postcontrol
for qq=length(jj);
    subplot(1,length(jj),qq);
    s3 = scatter_kde(y+(10+(rand(length(post_bird),1))*30),post_bird);
    s3.SizeData = 15;
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    xlabel('POST');
    [gg,bin]=hist(post_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
    hold on
    [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',8,'Color','r' ); % first peak
end

%% Density plot (all) - joy version %%

MM=[];
ff=[];

for wqw=1:41
    if (sum(data_matrix(:,wqw), 'omitnan')>0) == 0
        f=zeros(1,100);
        xi=zeros(1,100);
    else
        [f, Xi, u] = ksdensity(data_matrix(:,wqw));
    end
MM=[MM; Xi];
ff=[ff; f];
end
 
max(f);
x = linspace(0,9000,size(ff',1));

% figure(8)
% joyPlot(ff',x,0.0004,...
%     'StrokeColor','interp','StrokeWidth',3,'LineColor','none');
%     set(gca,'Color',[0.93 0.93 0.93]);
%     yticks([0.0004:0.0004:39])
% yticklabels({'POST',convertCharsToStrings(jj(end)),convertCharsToStrings(jj(end-1)),convertCharsToStrings(jj(end-2)),...
%     convertCharsToStrings(jj(end-3)),convertCharsToStrings(jj(end-4)),convertCharsToStrings(jj(end-5)),...
%     convertCharsToStrings(jj(end-6)),convertCharsToStrings(jj(end-7)),convertCharsToStrings(jj(end-8)),...
%     convertCharsToStrings(jj(end-9)),convertCharsToStrings(jj(end-10)),convertCharsToStrings(jj(end-11)),...
%     convertCharsToStrings(jj(end-12)),convertCharsToStrings(jj(end-13)),convertCharsToStrings(jj(end-14)),...
%     convertCharsToStrings(jj(end-15)),convertCharsToStrings(jj(end-16)),convertCharsToStrings(jj(end-17)),...
%     convertCharsToStrings(jj(end-18)),convertCharsToStrings(jj(end-19)),convertCharsToStrings(jj(end-20)),...
%      convertCharsToStrings(jj(end-21)),convertCharsToStrings(jj(end-22)),convertCharsToStrings(jj(end-23)),...
%     convertCharsToStrings(jj(end-24)),convertCharsToStrings(jj(end-25)),convertCharsToStrings(jj(end-26)),...
%     convertCharsToStrings(jj(end-27)),convertCharsToStrings(jj(end-28)),convertCharsToStrings(jj(end-29)),...
%     convertCharsToStrings(jj(end-30)),convertCharsToStrings(jj(end-31)),convertCharsToStrings(jj(end-32)),...
%     convertCharsToStrings(jj(end-33)),convertCharsToStrings(jj(end-34)),convertCharsToStrings(jj(end-35)),...
%     convertCharsToStrings(jj(end-36)),convertCharsToStrings(jj(end-37)),convertCharsToStrings(jj(end-38)),...
%     convertCharsToStrings(jj(end-39)),'PRE'})
%     ylabel('Bird freq (Hz)');
%     xlabel('Stimulus freq (Hz)');
%   % set(gca,'XGrid','on')
%     colormap((cool)); % Change the colormap, you know, for fun.
%     colorbar

%% Density plot (all) - 3D ribbon version %%
% figure(9)
% ribbon(ff')
%     zlabel('Bird freq (Hz)');
%     ylabel('Stimulus freq (Hz)');

%% Violin plot for quartiles %%

% meanPRE = mean(pre_bird,'omitnan');
% medianPRE = median(pre_bird,'omitnan');
% stdPRE = std(pre_bird,'omitnan');
% medianPRE+stdPRE;
% medianPRE-stdPRE;
% 
% figure
% violinplot(pre_bird);

%% Whistle characteristics %%
%% Duration of single whistles %%

post_bird=[post1_bird; post2_bird];
pre_bird=[pre1_bird; pre2_bird];

x=ones(length(pre_bird),1).*(1+(rand(length(pre_bird),1)-0.4)/5);
y=ones(length(post_bird),1).*(1+(rand(length(post_bird),1)-0.4)/5);
jj=[600:200:8600];

duration_PRE=[pre1_dur ;pre2_dur]/44100;

figure(10)
subplot(2,2,1)
violinplot(duration_PRE)
    ylabel('Whistle duration (s)');
    set(gca,'XTick',[]);
    xlim([0.5 1.5]);
    axis square
    
subplot(2,2,2)

histogram(duration_PRE,300,...
    'FaceColor', [80 160 50]/256, 'EdgeColor', [80 160 50]/256)
    ylabel('No. of whistles');
    xlabel('Whistle duration (s)');
    axis square

subplot(2,2,3)
scatter_kde(x,[duration_PRE]);
    hold on
    [f, Xi, u] = ksdensity(duration_PRE);
  % [a,b]=max(f)
    [a,b]=findpeaks(f);
    c=[a;b];
    cc=sortrows(c');
    freq1=Xi(cc(end,2));
    freq2=Xi(cc(end-1,2));
%    freq3=Xi(cc(end-2,2));
    yline(freq1,'-','LineWidth',8,'Color', 'r' ) % first peak
    yline(freq2,'-','LineWidth',8,'Color', 'm' ) % second peak
%    yline(freq3,'-','LineWidth',8,'Color', 'g' ) % third peak
    ylabel('Whistle duration (s)');
    set(gca,'XTick',[]);
    axis square

%% Duration of gaps between whistles %%

wgapsPRE = -[wgap1;wgap2] / 44100;
wgapsPRE = wgapsPRE(wgapsPRE > 0);
wgapsPRE = wgapsPRE(wgapsPRE < 0.5);
ggg = ones(length(wgapsPRE),1).*(1+(rand(length(wgapsPRE),1)-0.4) / 5);

figure(11)
subplot(2,2,1)
violinplot(wgapsPRE)
    axis square
    set(gca,'XTick',[]);
    xlim([0.5 1.5]);
    ylabel('Gap duration (s)');

subplot(2,2,2)

histogram(wgapsPRE,300,...
    'FaceColor', [80 160 50]/256, 'EdgeColor', [80 160 50]/256)
    ylabel('No. of whistles');
    xlabel('Gap duration (s)');
    axis square

subplot(2,2,3)
scatter_kde(ggg,[wgapsPRE]);
    hold on
    [f, Xi, u] = ksdensity(wgapsPRE);
  % [a,b]=max(f)
    [a,b]=findpeaks(f);
    c=[a;b];
    cc=sortrows(c');
    freq1=Xi(cc(end,2));
    freq2=Xi(cc(end-1,2));
  % freq3=Xi(cc(end-2,2));
    yline(freq1,'-','LineWidth',8,'Color', 'r'); % first peak
    yline(freq2,'-','LineWidth',8,'Color', 'm'); % second peak
  % yline(freq3,'-','LineWidth',8,'Color', 'k'); % third peak
    axis square
    ylabel('Gap duration (s)');
    set(gca,'XTick',[]);

%% Distribution of tempo %%

tempo = [tempo1; tempo2] / 44100;
tempo = tempo(tempo > 0);
tempo = tempo(tempo < 1.5);
ttt = ones(length(tempo),1).*(1+(rand(length(tempo),1)-0.4) / 5);

figure(12)
subplot(2,2,1)
violinplot(tempo)
    axis square
    ylabel('Tempo (s)');
    set(gca,'XTick',[]);
    xlim([0.5 1.5]);


subplot(2,2,2)
histogram(tempo,200)
    ylabel('No. of whistles');
    xlabel('Tempo (s)');

subplot(2,2,3)
scatter_kde(ttt,[tempo]);
    hold on
    [f, Xi, u] = ksdensity(tempo);
  % [a,b]=max(f)
    [a,b]=findpeaks(f);
    c=[a;b];
    cc=sortrows(c');
    freq1=Xi(cc(end,2));
    freq2=Xi(cc(end-1,2));
    freq3=Xi(cc(end-2,2));
    yline(freq1,'-','LineWidth',8,'Color', 'r' ); % first peak
    yline(freq2,'-','LineWidth',8,'Color', 'm' ); % second peak
    yline(freq3,'-','LineWidth',8,'Color', 'g' ); % third peak
    axis square
    ylabel('Tempo (s)');
    set(gca,'XTick',[]);


%% Density plot for freq 1000-1600 and PRE POST %%

post_bird=[post1_bird; post2_bird];
pre_bird=[pre1_bird; pre2_bird];

x=ones(length(pre_bird),1).*(1+(rand(length(pre_bird),1)-0.4)/5);
y=ones(length(post_bird),1).*(1+(rand(length(post_bird),1)-0.4)/5);
www=[800:200:1600];

data_matrix=nan(10000);
data_matrix(1:length(pre_bird),1)=pre_bird;
data_matrix(1:length(post_bird),length(www)+1)=post_bird;

for w=2:length (www)
    ids=find(pb_stim==www(w));
    data_matrix(1:length(ids),w)=pb_bird(ids);
end

whistle_under_pref=[];
stim_under_pref=[];
peak_low=[];
figure(13)
for qq=1
    subplot(1,length(www)+1,qq)
    scatter_kde(x+(10+(rand(length(pre_bird),1))*30),pre_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('PRE');
    hold on
    [gg,bin]=hist(pre_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' );
    [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' );
end

for qq=2:length(www);
    subplot(1,length(www)+1,qq);
    ids=find(pb_stim==www(qq));
    pb_bird(ids);
     zzz=pb_bird(ids);
    www(qq);
    [p,h]=ranksum(zzz,www(qq));
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids)+(10+(rand(length(ids),1))*30),pb_bird(ids));
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(www(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
       %yline(median(pb_bird(ids),'omitnan'),'-','Color', 'r', 'LineWidth',8);
       %yline(mean(pb_bird(ids),'omitnan'),'-','Color', 'm', 'LineWidth',8);
        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
       %yline(pb_max,'-','LineWidth',8,'Color', 'r' )
        hold on
        [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
       %[a,b]=max(f)
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
       % freq2=Xi(cc(end-1,2));
        yline(freq1,'-','LineWidth',8,'Color', 'r' );
      % yline(freq2,'-','LineWidth',8,'Color', 'm' );
      peak_low=[peak_low;freq1];
    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(www(qq)));
        hold on
    end
    hold on
    
   
    
    stim_under_pref=[stim_under_pref; pb_stim(ids)];
    whistle_under_pref=[whistle_under_pref; pb_bird(ids)];

end

for qq=length(www)+1;
    subplot(1,length(www)+1,qq);
    scatter_kde(y+(10+(rand(length(post_bird),1))*30),post_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('POST');
    [gg,bin]=hist(post_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
    hold on
    [f, Xi, u] = ksdensity(post_bird,'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' )
end

%%  Density plot for freq under pref range only!  and PRE POST %%

www=[800:200:1600];

whistle_under_pref=[];
stim_under_pref=[];
peak_low_only=[];
figure(1343)
for qq=1
    subplot(1,length(www)+1,qq)
    scatter_kde(x+(10+(rand(length(pre_bird),1))*30),pre_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('PRE');
    hold on
    [gg,bin]=hist(pre_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' );
    [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' );
end

for qq=2:length(www);
   
   subplot(1,length(www)+1,qq);
    ids=find(pb_stim==www(qq));
    pb_bird(ids);
      ddd=pb_bird(ids);
     ix=find(ddd>3700);
     l=ddd(ix);
    ixx=find(ddd<1700);
     ll=ddd(ixx);
    lll=[ l;ll];
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids(1:length(lll)))+(10+(rand(length(ids(1:length(lll))),1))*30),lll);
        ylim([0 9000]);
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(www(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
       %yline(median(pb_bird(ids),'omitnan'),'-','Color', 'r', 'LineWidth',8);
       %yline(mean(pb_bird(ids),'omitnan'),'-','Color', 'm', 'LineWidth',8);
        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
       %yline(pb_max,'-','LineWidth',8,'Color', 'r' )
        hold on
        [f, Xi, u] = ksdensity(lll,'Support','positive'); %data_matrix(:,qq));
       %[a,b]=max(f)
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        %freq2=Xi(cc(end-1,2));
        yline(freq1,'-','LineWidth',8,'Color', 'r' );
      % yline(freq2,'-','LineWidth',8,'Color', 'm' );
      peak_low_only=[peak_low_only;freq1];
    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(www(qq)));
        hold on
    end
    hold on
    
%      figure
% subplot(1,2,1)
%     h=scatterhist(pb_stim(ids(1:length(lll)))+(10+(rand(length(ids(1:length(lll))),1))*30),lll,'Location','SouthEast',...
   % 'Direction','out')
%     delete(h(2))
%     yline(1200)
%     ylim([0 9000])
%   
 

     stim_under_pref=[stim_under_pref; pb_stim(ids(1:length(lll)))];
    whistle_under_pref=[whistle_under_pref; lll];
   

end

for qq=length(www)+1;
    subplot(1,length(www)+1,qq);
    scatter_kde(y+(10+(rand(length(post_bird),1))*30),post_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('POST');
    [gg,bin]=hist(post_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
    hold on
    [f, Xi, u] = ksdensity(post_bird,'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' )
end


%%
tpb = table(whistle_under_pref, stim_under_pref);
mpb = fitlm(tpb,'whistle_under_pref ~ stim_under_pref')

figure(14)
ppb =plot(mpb,'color','k','Marker','o','MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','k','LineWidth',3);
    ylabel('Bird freq (Hz)')
    xlabel('Stimulus freq (Hz)')
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
    title off
    xlim([500 2000]);
    ylim([0 9000]);
    title 'Freq. bird VS stimuli (1000-1600 Hz)';
          
%% Density plot for freq 1800-3400 and PRE POST %%

post_bird=[post1_bird; post2_bird];
pre_bird=[pre1_bird; pre2_bird];

x=ones(length(pre_bird),1).*(1+(rand(length(pre_bird),1)-0.4)/5);
y=ones(length(post_bird),1).*(1+(rand(length(post_bird),1)-0.4)/5);
ww=[1600:200:3600];

data_matrix=nan(10000);
data_matrix(1:length(pre_bird),1)=pre_bird;
data_matrix(1:length(post_bird),41)=post_bird;

for w=2:length (ww)
    ids=find(pb_stim==ww(w));
    data_matrix(1:length(ids),w)=pb_bird(ids);
end

whistle_under_pref=[];
stim_under_pref=[];

figure(15)
for qq=1
    subplot(1,length(ww)+1,qq)
    scatter_kde(x+(10+(rand(length(pre_bird),1))*30),pre_bird);
    xlim([-10 50]);
    ylim([1000 9000]);
    set(gca,'xtick',[]);
    xlabel('PRE');
    hold on
    [gg,bin]=hist(pre_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' );
    [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' );
    set(gca,'TickDir','out')
set(gca, 'yScale', 'log')

end

peak_inside=[];
for qq=2:length(ww);
    subplot(1,length(ww)+1,qq);
    ids=find(pb_stim==ww(qq));
    pb_bird(ids);
    zzz=pb_bird(ids);
    ww(qq);
    %[p,h]=ranksum(zzz,pre_bird);
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids)+(10+(rand(length(ids),1))*30),pb_bird(ids));
        ylim([1000 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(ww(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
       %yline(median(pb_bird(ids),'omitnan'),'-','Color', 'r', 'LineWidth',8);
       %yline(mean(pb_bird(ids),'omitnan'),'-','Color', 'm', 'LineWidth',8);
        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
       %yline(pb_max,'-','LineWidth',8,'Color', 'r' )
        hold on
        [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
       %[a,b]=max(f)
        [a,b]=findpeaks(f);
        c=[a;b];
         cc=sortrows(c');
        freq1=Xi(cc(end,2));
        %freq2=Xi(cc(end-1,2));
        yline(freq1,'-','LineWidth',8,'Color', 'r' );
       %yline(freq2,'-','LineWidth',8,'Color', 'm' );
       peak_inside=[peak_inside;freq1];
    else 
        ylim([1000 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(ww(qq)));
        hold on
    end
    hold on
    set(gca, 'yScale', 'log')

end

for qq=length(ww)+1;
    subplot(1,length(ww)+1,qq);
    scatter_kde(y+(10+(rand(length(post_bird),1))*30),post_bird);
    xlim([-10 50]);
    ylim([1000 9000]);
    set(gca,'xtick',[]);
    xlabel('POST');
    [gg,bin]=hist(post_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
    hold on
    [f, Xi, u] = ksdensity(post_bird,'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' )
        set(gca,'TickDir','out')
set(gca,'YAxisLocation', 'right')
set(gca, 'yScale', 'log')

end
%% inside preferred - only preferred!

 stim_pref=[];
    whistle_pref=[];

figure(1345)
for qq=1
    subplot(1,length(ww)+1,qq)
    scatter_kde(x+(10+(rand(length(pre_bird),1))*30),pre_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('PRE');
    hold on
    [gg,bin]=hist(pre_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' );
    [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
    [a,b]=max(f);
    freqPRE=Xi(b);
    yline(freqPRE,'-','LineWidth',8,'Color', 'r' );
end

peak_inside_only=[];
for qq=2:length(ww);
    subplot(1,length(ww)+1,qq);
    ids=find(pb_stim==ww(qq));
    pb_bird(ids);
    ddd=pb_bird(ids);
    ix=find(ddd<3700);
    l=ddd(ix);
    ixx=find(ddd>1600);
     ll=ddd(ixx);
    lll=intersect(l ,ll);
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids(1:length(lll)))+(10+(rand(length(ids(1:length(lll))),1))*30),lll);
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(ww(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
       %yline(median(pb_bird(ids),'omitnan'),'-','Color', 'r', 'LineWidth',8);
       %yline(mean(pb_bird(ids),'omitnan'),'-','Color', 'm', 'LineWidth',8);
        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
       %yline(pb_max,'-','LineWidth',8,'Color', 'r' )
        hold on
        [f, Xi, u] = ksdensity(lll,'Support','positive');%data_matrix(:,qq));
       %[a,b]=max(f)
        [a,b]=findpeaks(f);
        c=[a;b];
         cc=sortrows(c');
        freq1=Xi(cc(end,2));
        %freq2=Xi(cc(end-1,2));
        yline(freq1,'-','LineWidth',8,'Color', 'r' );
       %yline(freq2,'-','LineWidth',8,'Color', 'm' );
       peak_inside_only=[peak_inside_only;freq1];
    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(ww(qq)));
        hold on
    end
    hold on
      stim_pref=[stim_pref; pb_stim(ids(1:length(lll)))];
    whistle_pref=[whistle_pref; lll];
end

for qq=length(ww)+1;
    subplot(1,length(ww)+1,qq);
    scatter_kde(y+(10+(rand(length(post_bird),1))*30),post_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('POST');
    [gg,bin]=hist(post_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
    hold on
    [f, Xi, u] = ksdensity(post_bird,'Support','positive');
    [a,b]=max(f);
    freq=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' )
end


%%
tpb = table(whistle_pref, stim_pref);
mpb = fitlm(tpb,'whistle_pref ~ stim_pref')

figure(5745)
hold on
ppb =plot(mpb,'color','m','Marker','o','MarkerFaceColor', 'm', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','k','LineWidth',3);
    ylabel('Bird freq (Hz)')
    xlabel('Stimulus freq (Hz)')
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
    title off
    xlim([1400 3800]);
    ylim([0 9000]);
    title 'Freq. bird VS stimuli (1800-3400 Hz)';

 %% Density plot for freq 3800-8000 and PRE POST %%
post_bird=[post1_bird; post2_bird];
pre_bird=[pre1_bird; pre2_bird];

x=ones(length(pre_bird),1).*(1+(rand(length(pre_bird),1)-0.4)/5);
y=ones(length(post_bird),1).*(1+(rand(length(post_bird),1)-0.4)/5);
wwx=[3600:200:9000];

data_matrix=nan(10000);
data_matrix(1:length(pre_bird),1)=pre_bird;
data_matrix(1:length(post_bird),length(wwx)+1)=post_bird;

for w=2:length (wwx)
    ids=find(pb_stim==wwx(w));
    data_matrix(1:length(ids),w)=pb_bird(ids);
 
end

whistle_over_pref=[];
stim_over_pref=[];
peak_high=[];
figure(17)
for qq=1
    subplot(1,length(wwx)+1,qq)
    scatter_kde(x+(10+(rand(length(pre_bird),1))*30),pre_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('PRE');
    hold on
    [gg,bin]=hist(pre_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' );
    [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
    [a,b]=max(f);
    freq_pre=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' );
end
second_peak_high=[];
for qq=2:length(wwx);
    subplot(1,length(wwx)+1,qq);
    ids=find(pb_stim==wwx(qq));
    pb_bird(ids);
    l=pb_bird(pb_bird(ids)>3700);
    ll=pb_bird(pb_bird(ids)<1700);
    lll=[ll; l];
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids)+(10+(rand(length(ids),1))*30),pb_bird(ids));
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwx(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
       %yline(median(pb_bird(ids),'omitnan'),'-','Color', 'r', 'LineWidth',8);
       %yline(mean(pb_bird(ids),'omitnan'),'-','Color', 'm', 'LineWidth',8);
        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
       %yline(pb_max,'-','LineWidth',8,'Color', 'r' )
        hold on
        [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
       %[a,b]=max(f)
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        %freq2=Xi(cc(end-1,2));
       %second_peak_high=[second_peak_high, freq2];
        yline(freq1,'-','LineWidth',8,'Color', 'r' );
       % yline(freq2,'-','LineWidth',8,'Color', 'm' );
               peak_high=[peak_high;freq1];

    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwx(qq)));
        hold on
    end
    hold on
    stim_over_pref=[stim_over_pref; pb_stim(ids)];
    whistle_over_pref=[whistle_over_pref; pb_bird(ids)];
end

for qq=length(wwx)+1;
    subplot(1,length(wwx)+1,qq);
    scatter_kde(y+(10+(rand(length(post_bird),1))*30),post_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('POST');
    [gg,bin]=hist(post_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
    % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
    hold on
    [f, Xi, u] = ksdensity(post_bird,'Support','positive');
    [a,b]=max(f);
    freq_post=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' )
end

%%
post_bird=[post1_bird; post2_bird];
pre_bird=[pre1_bird; pre2_bird];

x=ones(length(pre_bird),1).*(1+(rand(length(pre_bird),1)-0.4)/5);
y=ones(length(post_bird),1).*(1+(rand(length(post_bird),1)-0.4)/5);
wwx=[3600:200:9000];

data_matrix=nan(10000);
data_matrix(1:length(pre_bird),1)=pre_bird;
data_matrix(1:length(post_bird),length(wwx)+1)=post_bird;

for w=2:length (wwx)
    ids=find(pb_stim==wwx(w));
    data_matrix(1:length(ids),w)=pb_bird(ids);
 
end

whistle_over_pref=[];
stim_over_pref=[];
peak_high_only=[];
figure(117)
for qq=1
    subplot(1,length(wwx)+1,qq)
    scatter_kde(x+(10+(rand(length(pre_bird),1))*30),pre_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('PRE');
    hold on
    [gg,bin]=hist(pre_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
  % yline(pb_max,'-','LineWidth',8,'Color', 'r' );
    [f, Xi, u] = ksdensity(data_matrix(:,qq),'Support','positive');
    [a,b]=max(f);
    freq_pre=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' );
end

for qq=2:length(wwx);
    subplot(1,length(wwx)+1,qq);
    ids=find(pb_stim==wwx(qq));
    ddd=pb_bird(ids);
    ix=find(ddd>3652);
     l=ddd(ix);
     ixx=find(ddd<1700);
      ll=ddd(ixx);
    lll=[l];
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids(1:length(lll)))+(10+(rand(length(ids(1:length(lll))),1))*30),lll);
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwx(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
       %yline(median(pb_bird(ids),'omitnan'),'-','Color', 'r', 'LineWidth',8);
       %yline(mean(pb_bird(ids),'omitnan'),'-','Color', 'm', 'LineWidth',8);
        [gg,bin]=hist(lll,25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
       %yline(pb_max,'-','LineWidth',8,'Color', 'r' )
        hold on
        [f, Xi, u] = ksdensity(lll,'Support','positive');%data_matrix(:,qq));
       %[a,b]=max(f)
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        %freq2=Xi(cc(end-1,2));
              %yline(freq2,'-','LineWidth',8,'Color', 'm' );

        yline(freq1,'-','LineWidth',8,'Color', 'r' );
        %yline(freq2,'-','LineWidth',8,'Color', 'm' );
               peak_high_only=[peak_high_only;freq1];

    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwx(qq)));
        hold on
    end
    hold on
    stim_over_pref=[stim_over_pref; pb_stim(ids(1:length(lll)))];
    whistle_over_pref=[whistle_over_pref; lll];
end



for qq=length(wwx)+1;
    subplot(1,length(wwx)+1,qq);
    scatter_kde(y+(10+(rand(length(post_bird),1))*30),post_bird);
    xlim([-10 50]);
    ylim([0 9000]);
    set(gca,'xtick',[]);
    xlabel('POST');
    [gg,bin]=hist(post_bird);
    [r_max,b]=max(gg);
    pb_max=bin(b);
    % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
    hold on
    [f, Xi, u] = ksdensity(post_bird,'Support','positive');
    [a,b]=max(f);
    freq_post=Xi(b);
    yline(freq,'-','LineWidth',8,'Color', 'r' )
end


%%
tpb = table(whistle_over_pref, stim_over_pref);
mpb = fitlm(tpb,'whistle_over_pref ~ stim_over_pref')
rr=[];

for vv=1:29;
     if (sum(data_matrix(:,vv), 'omitnan')>0)==0;
        r=zeros(1,100);
        xi=zeros(1,100);
     else
        [r, Xi, u] = ksdensity(data_matrix(:,vv));
     end
    MM=[MM; Xi];
    rr=[rr; r];
end
 
figure(18)
ppb = plot(mpb,'color','k','Marker','o','MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','k','LineWidth',3);
    ylabel('Bird freq (Hz)')
    xlabel('Stimulus freq (Hz)')
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
    title off
    xlim([ 3400 9000]);
    ylim([0 9000]);
    title 'Frequency stimuli VS bird (3800-9000 Hz)'
     %%
% figure(19)
% subplot(1,3,1)
% ppb = plot(mpb,'color','k','Marker','o','MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
%     hold on
%     set(ppb(2), 'Color','k','LineWidth',3);
%     ylabel('Bird freq (Hz)')
%     xlabel('Stimulus freq (Hz)')
%     box off
%     set(gca,'linewidth',1,'FontSize', 14)
%     legend off
%     title off
%     xlim([ 3400 9000]);
%     ylim([0 9000]);
%     title 'Frequency stimuli VS bird (3800-9000 Hz)';
%     axis square
% 
% subplot(1,3,2)
% ribbon(rr')
%     zlabel('Bird freq (Hz)');
%     ylabel('Stimulus freq (Hz)');
%     axis square
% 
% subplot(1,3,3)
% xxxx = linspace(0,9000,size(rr',1));
% joyPlot(rr',xxxx,0.0004,...
%     'StrokeColor','interp','StrokeWidth',3,'LineColor','none')
%     set(gca,'Color',[0.93 0.93 0.93])
%     yticklabels({'POST',convertCharsToStrings(jj(end)),convertCharsToStrings(jj(end-1)),convertCharsToStrings(jj(end-2)),...
%     convertCharsToStrings(jj(end-3)),convertCharsToStrings(jj(end-4)),convertCharsToStrings(jj(end-5)),...
%     convertCharsToStrings(jj(end-6)),convertCharsToStrings(jj(end-7)),convertCharsToStrings(jj(end-8)),...
%     convertCharsToStrings(jj(end-9)),convertCharsToStrings(jj(end-10)),convertCharsToStrings(jj(end-11)),...
%     convertCharsToStrings(jj(end-12)),convertCharsToStrings(jj(end-13)),convertCharsToStrings(jj(end-14)),...
%     convertCharsToStrings(jj(end-15)),convertCharsToStrings(jj(end-16)),convertCharsToStrings(jj(end-17)),...
%     convertCharsToStrings(jj(end-18)),convertCharsToStrings(jj(end-19)),convertCharsToStrings(jj(end-20)),...
%      convertCharsToStrings(jj(end-21)),convertCharsToStrings(jj(end-22)),convertCharsToStrings(jj(end-23)),...
%     convertCharsToStrings(jj(end-24)),convertCharsToStrings(jj(end-25)),convertCharsToStrings(jj(end-26)),...
%     convertCharsToStrings(jj(end-27)),convertCharsToStrings(jj(end-28)),convertCharsToStrings(jj(end-29)),...
%     convertCharsToStrings(jj(end-30)),convertCharsToStrings(jj(end-31)),convertCharsToStrings(jj(end-32)),...
%     convertCharsToStrings(jj(end-33)),convertCharsToStrings(jj(end-34)),convertCharsToStrings(jj(end-35)),...
%     convertCharsToStrings(jj(end-36)),convertCharsToStrings(jj(end-37)),convertCharsToStrings(jj(end-38)),...
%     convertCharsToStrings(jj(end-39)),'PRE'})
%     ylabel('Bird freq (Hz)');
%     xlabel('Stimulus freq (Hz)');
%     colormap(flipud(copper));
%     colorbar
%     axis square
%%
figure

for qq=1:5;
    subplot(1,length(jj(2:end-2)),qq);
     ids=find(pb_stim==www(qq));
    pb_bird(ids);
      ddd=pb_bird(ids);
     ix=find(ddd>3700);
     l=ddd(ix);
    ixx=find(ddd<1700);
     ll=ddd(ixx);
    lll=[ l;ll];
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids(1:length(lll)))+(10+(rand(length(ids(1:length(lll))),1))*30),lll);
        ylim([1000 9000]);
       
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(www(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
        hold on
        [f, Xi, u] = ksdensity(lll,'Support','positive');
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        yline(freq1,'-','LineWidth',8,'Color', 'r' );
    else 
        ylim([1000 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(www(qq)));
        hold on
    end
    hold on
         set(gca, 'yScale', 'log')

end

wwww=jj(2:end-2);

for qq=6:15;
      subplot(1,length(jj(2:end-2)),qq);
    ids=find(pb_stim==wwww(qq));
    pb_bird(ids);
    zzz=pb_bird(ids);
    wwww(qq);
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids)+(10+(rand(length(ids),1))*30),pb_bird(ids));
        ylim([1000 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwww(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
        hold on
        [f, Xi, u] = ksdensity(zzz,'Support','positive');
        [a,b]=findpeaks(f);
        c=[a;b];
         cc=sortrows(c');
        freq1=Xi(cc(end,2));
        yline(freq1,'-','LineWidth',8,'Color', 'r' );
    else 
        ylim([1000 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwww(qq)));
        hold on
    end
    hold on
     set(gca, 'yScale', 'log')
 
end

for qq=16:length(jj(2:end-2));
    subplot(1,length(jj(2:end-2)),qq);
    ids=find(pb_stim==wwww(qq));
    ddd=pb_bird(ids);
    ix=find(ddd>3700);
    l=ddd(ix);
    ixx=find(ddd<1700);
     ll=ddd(ixx);
    lll=[l; ll];
    if (sum(pb_stim(ids))>0)==1;
        scatter_kde(pb_stim(ids(1:length(lll)))+(10+(rand(length(ids(1:length(lll))),1))*30),lll);
        ylim([1000 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwww(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',8);
        [gg,bin]=hist(lll,25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
        hold on
        [f, Xi, u] = ksdensity(lll,'Support','positive');
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        yline(freq1,'-','LineWidth',8,'Color', 'r' );
    else 
        ylim([1000 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwww(qq)));
        hold on
    end
    hold on
         set(gca, 'yScale', 'log')

end


%%

outside_bird=[peak_low; peak_high];
outside_pb=[www(2:end) [3800,4000,4400,5000,6000,7000,8000]];

high=[3800,4000,4400,5000,6000,7000,8000];

inside_bird=peak_inside;
inside_pb=(ww(2:end))';

[p,h]= signrank(outside_pb, outside_bird)
[p,h] =signrank(inside_pb, inside_bird)
[p,h] =signrank(high, second_peak_high)

figure
subplot(1,3,1)
plot([2*ones(1,length(outside_pb))], (outside_bird),...
    'o','MarkerFaceColor', 'r', 'MarkerEdgeColor','r', 'MarkerSize', 7)
hold on
plot([1*ones(1,length(outside_bird))], (outside_pb),...
    'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 7)
hold on
plot([1*ones(1,length(outside_bird));2* ones(1,length(outside_bird))],[outside_pb',outside_bird]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on
xlim([0 3])
ylim([0 9000])
ylabel('% of whistle songs')
box off
xticks([1 2 3 4])
xticklabels({'Pairs','PRE','PB','POST'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on

subplot(1,3,2)
plot([2*ones(1,length(inside_pb))], (inside_bird),...
    'o','MarkerFaceColor', 'r', 'MarkerEdgeColor','r', 'MarkerSize', 7)
hold on
plot([1*ones(1,length(inside_bird))], (inside_pb),...
    'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 7)
hold on
plot([1*ones(1,length(inside_bird));2* ones(1,length(inside_bird))],[inside_pb,inside_bird]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on
hold on
xlim([0 3])
ylim([0 9000])
ylabel('% of whistle songs')
box off
xticks([1 2 3 4])
xticklabels({'Pairs','PRE','PB','POST'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on


% subplot(1,3,3)
% plot([1*ones(1,length(high))], (high),...
%     'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 7)
% hold on
% plot([2*ones(1,length(high))], (second_peak_high),...
%     'o','MarkerFaceColor', 'm', 'MarkerEdgeColor','m', 'MarkerSize', 7)
% hold on
% plot([1*ones(1,length(high));2* ones(1,length(high))],[high' second_peak_high']',...
%     '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
% hold on
% hold on
% xlim([0 3])
% ylim([0 9000])
% ylabel('% of whistle songs')
% box off
% xticks([1 2 3 4])
% xticklabels({'Pairs','PRE','PB','POST'})
% %legend(sprintf('p = %0.3f',p_anti))
% %title('Percentage of whistle songs')
% axis square
% set(gca,'linewidth',1,'FontSize', 14)
% hold on
%%
outside_only_bird=[peak_low_only; peak_high_only];
outside_pb=[www(2:end) [3800,4000,4400,5000,6000,7000,8000]];
low=www(2:end)
high=[3800,4000,4400,5000,6000,7000,8000];

inside_only_bird=peak_inside_only;
inside_pb=(ww(2:end))';


figure
subplot(1,3,1)
plot([2*ones(1,length(low))], (peak_low_only),...
    'o','MarkerFaceColor', 'r', 'MarkerEdgeColor','r', 'MarkerSize', 7)
hold on
plot([1*ones(1,length(peak_low_only))], (low),...
    'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 7)
hold on
plot([1*ones(1,length(low));2* ones(1,length(peak_low_only))],[low' peak_low_only]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on
xlim([0 3])
ylim([0 9000])
ylabel('% of whistle songs')
box off
xticks([1 2 3 4])
xticklabels({'Pairs','PRE','PB','POST'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on

[p,h]= signrank(low' ,peak_low_only)

subplot(1,3,2)
plot([2*ones(1,length(inside_pb))], (inside_bird),...
    'o','MarkerFaceColor', 'r', 'MarkerEdgeColor','r', 'MarkerSize', 7)
hold on
plot([1*ones(1,length(inside_bird))], (inside_pb),...
    'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 7)
hold on
plot([1*ones(1,length(inside_bird));2* ones(1,length(inside_bird))],[inside_pb,inside_bird]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on
hold on
xlim([0 3])
ylim([0 9000])
ylabel('% of whistle songs')
box off
xticks([1 2 3 4])
xticklabels({'Pairs','PRE','PB','POST'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on

[p,h] =signrank(inside_pb,inside_only_bird)

subplot(1,3,3)
plot([1*ones(1,length(high))], (high),...
    'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 7)
hold on
plot([2*ones(1,length(high))], (peak_high_only),...
    'o','MarkerFaceColor', 'r', 'MarkerEdgeColor','r', 'MarkerSize', 7)
hold on
plot([1*ones(1,length(high));2* ones(1,length(high))],[high' peak_high_only]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on
hold on
xlim([0 3])
ylim([0 9000])
ylabel('% of whistle songs')
box off
xticks([1 2 3 4])
xticklabels({'Pairs','PRE','PB','POST'})
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
set(gca,'linewidth',1,'FontSize', 14)
hold on

[p,h] =signrank(high', peak_high_only)

[p,h] =signrank([low'; high'], [peak_low_only;peak_high_only])


%%
figure
subplot(1,2,1)
plot([2*ones(1,length(outside_bird))], (outside_pb),...
    'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([1*ones(1,length(outside_pb))], (outside_bird),...
    'o','MarkerFaceColor', [255 15 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([3*ones(1,length(inside_pb))], (inside_bird),...
    'o','MarkerFaceColor', [255 15 65]/256, 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
plot([2*ones(1,length(inside_bird))], (inside_pb),...
    'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[215 105 0]/256, 'MarkerSize', 7)
hold on
hold on
plot([2*ones(1,length(outside_bird));1* ones(1,length(outside_bird))],[outside_pb',outside_bird]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on
plot([2*ones(1,length(inside_bird));3* ones(1,length(inside_bird))],[inside_pb,inside_bird]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
xlim([0 6])
ylim([0 9000])

axis square

ylabel('pitch (Hz)')

plot([5*ones(1,length(low))], (peak_low_only),...
    'o','MarkerFaceColor', 'r', 'MarkerEdgeColor','r', 'MarkerSize', 7)
hold on
hold on
plot([2*ones(1,length(low));5* ones(1,length(peak_low_only))],[low' peak_low_only]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)

plot([5*ones(1,length(high))], (peak_high_only),...
    'o','MarkerFaceColor', 'r', 'MarkerEdgeColor','r', 'MarkerSize', 7)
hold on
plot([2*ones(1,length(high));5* ones(1,length(high))],[high' peak_high_only]',...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)

xticks([1 2 3,5])
xticklabels({'ouside preferred range','PB','inside preferred range', 'ouside and without preferred range'})
box off
set(gca,'TickDir','out')

axis square
%%
delta_all=peak_all_bird - peak_all_pb;

delta_outside=outside_bird - outside_pb';

delta_inside=inside_bird - inside_pb;
delta_inside_only=peak_inside_only - inside_pb;


delta_inside_pre=inside_bird - freqPRE;

delta_low=peak_low_only - low';

delta_low_pre=peak_low_only - freqPRE;

delta_high_pre=peak_high_only - freqPRE;

delta_high=peak_high_only - high';

delta=[delta_low;delta_inside  ;delta_high];
delta_pre=[delta_low_pre;delta_inside_pre;  delta_high_pre];

delta_all=[delta_outside(1:4);delta_inside;  delta_outside(5:end)];

%%
pbs=[low'; inside_pb; high'];
figure
subplot(2,2,1)
for qww=1:length(delta)
   e= [rand rand rand];
stem(pbs(qww),delta(qww),'filled','LineWidth',2,'LineStyle','-.','Color',[e],...
     'MarkerFaceColor',[e], 'MarkerSize', 8)
 hold on
end
xlim([0 9000])
ylim([-1000 1000])
xline(1700,'--')
xline(3700,'--')
box off
set(gca,'TickDir','out')
yticks([-1000 -500 0 500 1000])

%set(gca, 'XScale', 'log')

subplot(2,2,3)
for qww=1:length(delta)
   e= [rand rand rand];
stem(pbs(qww),delta(qww),'filled','LineWidth',2,'LineStyle','-.','Color',[e],...
     'MarkerFaceColor',[e], 'MarkerSize', 8)
 hold on
end
xlim([0 9000])
ylim([-5500 -3500])
xline(1700,'--')
xline(3700,'--')
box off

set(gca,'TickDir','out')
%set(gca, 'XScale', 'log')

% subplot(1,3,2)
% CategoricalScatterplot(delta,ones(length(delta),1))
% ylim([-6000 2000])
% xlim([0 2])
 subplot(2,2,2)
violinplot(delta)
ylim([-1000 1000])
xlim([0.5 1.5])
box off
set(gca,'TickDir','out')
yline(0, '-')
axis square
yticks([-1000 -500 0 500 1000])

subplot(2,2,4)
violinplot(delta)
ylim([-5500 -3500])
xlim([0.5 1.5])
box off
set(gca,'TickDir','out')
yline(0, '-')
axis square

[p,h]=ranksum(delta_outside,delta_inside)
std(delta)
%%
figure
raincloud_plot((delta_inside))
hold on
raincloud_plot((delta_outside))
xline(-0)
xlim([-8000 8000])
ylim([-0.0025 0.0025])


group1=ones(length(delta_inside),1);
group2=ones(length(delta_outside),1)*2;
groupp=[group1; group2];
anova1(([(delta_inside);(delta_outside)]),groupp)

group1=ones(length(delta_inside),1);
group2=ones(length(delta_inside_pre),1)*2;
groupp=[group1; group2];
anova1(([delta_inside;delta_inside_pre]),groupp)

group1=ones(length(delta_inside),1);
group2=ones(length(delta_outside),1)*2;
groupp=[group1; group2];

figure
violinplot(([delta_inside;delta_outside]),groupp)
xlim([0 3])
ylim([-6000 6000])
axis square
set(gca,'linewidth',1,'FontSize', 14)

figure
CategoricalScatterplot(([delta_inside;delta_outside]), groupp)
xlim([0 3])
ylim([-6000 6000])
axis square
set(gca,'linewidth',1,'FontSize', 14)

d=[delta_low; delta_high];
group1=ones(length(d),1);
group2=ones(length(delta_outside),1)*2;
groupp=[group1; group2];

figure
CategoricalScatterplot(([d;delta_outside]), groupp)
anova1(([d(1:end);delta_outside(1:end)]), groupp)

group1=ones(length(d),1);
group2=ones(length(delta_inside),1)*2;
groupp=[group1; group2];

figure
CategoricalScatterplot(([d;delta_inside]), groupp)
anova1(([abs(d(1:end));abs(delta_inside(1:end))]), groupp)


group1=ones(length(delta_inside),1);
group2=ones(length(delta_outside),1)*2;
group3=ones(length(d),1)*3;
grouppp=[group1; group2; group3];

figure
CategoricalScatterplot(abs([delta_inside;delta_outside;d]), grouppp)
anova1(([delta_inside;delta_outside]), [group1;group2])
anova1(([delta_inside;d]), [group1;group3])
anova1(abs([delta_outside;d]), [group2;group3])

median(d)


anova1(([delta_inside;delta_outside;d]), [group1;group2;group3])


figure
raincloud_plot([delta_inside;d])
xline(0,'--')
signrank([delta_inside;d],0)

[p,h]=ranksum(abs(delta_inside),abs(delta_outside))
[p,h]=ranksum(abs(delta_inside),abs(d))
[p,h]=ranksum(abs(delta_outside),abs(d))

p=ranksum(inside_bird, outside_bird)

median(inside_bird)
median(outside_bird)
%%

pbs=[low'; inside_pb; high'];
figure
subplot(2,1,1)
for qww=1:length(delta)
   e= [rand rand rand];
stem(pbs(qww),(delta_all(qww)),'filled','LineWidth',2,'LineStyle','-','Color',[e],...
     'MarkerFaceColor',[e], 'MarkerSize', 8)
 hold on
end
xlim([0 9000])
 ylim([-6000 3000])
%ylim([0 6000])
xline(1700,'--')
xline(3700,'--')
box off
set(gca,'TickDir','out')
% set(gca, 'yScale', 'log')
axis square
subplot(2,1,2)
for qww=1:length(delta)
   e= [rand rand rand];
stem(pbs(qww),(delta(qww)),'filled','LineWidth',2,'LineStyle','-.','Color',[e],...
     'MarkerFaceColor',[e], 'MarkerSize', 8)
 hold on
end
xlim([0 9000])
 ylim([-6000 3000])
%ylim([0 6000])
xline(1700,'--')
xline(3700,'--')
box off
set(gca,'TickDir','out')
% set(gca, 'yScale', 'log')

% subplot(1,3,2)
% CategoricalScatterplot(delta,ones(length(delta),1))
% ylim([-6000 2000])
% xlim([0 2])
 subplot(2,2,2)
CategoricalScatterplot(([delta_inside;delta_outside;[delta_inside;delta_outside ]]),[ones(length(delta_inside),1);2*ones(length(delta_outside),1)...
   ;3*ones(length([delta_inside;delta_outside ]),1)])
% CategoricalScatterplot([delta_inside;delta_outside],[ones(length(delta_inside),1);2*ones(length(delta_outside),1)])
ylim([-6000 3000])
xlim([0.5 2.5])
box off
set(gca,'TickDir','out')
yline(0, '-')
axis square

subplot(2,2,4)
 CategoricalScatterplot(([delta_inside;[delta_low ;delta_high];delta]),[ones(length(delta_inside),1);[2*ones(length(delta_low),1);2*ones(length(delta_high),1)];...
   3*ones(length(delta),1)])
% CategoricalScatterplot([delta_inside;[delta_low ;delta_high]],[ones(length(delta_inside),1);...
%     [2*ones(length(delta_low),1);2*ones(length(delta_high),1)]])
% 
ylim([-6000 6000])
xlim([0.5 2.5])
box off
set(gca,'TickDir','out')
yline(0, '-')
axis square

[p,h]=ranksum(delta_inside ,delta_outside(1:4))

[h,p]=ttest([delta_outside])

p=ranksum(delta_inside, 0)



[p,h]=ttest([delta_inside])


mean(([delta_inside]))

std(([delta_outside]))/100

median(abs([delta_low ;delta_high]))

[h,p]=lillietest(delta_outside)

%%

del=((delta_all))
 tpb = table(pbs, (del));
mpb = fitlm(tpb,' del~pbs ')
figure
ppb = plot(mpb,'color','k','Marker','none','MarkerSize',4,'MarkerFaceColor', 'none', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','k','LineWidth',3);
    legend off


 [values, centers] = hist3([   pbs  delta_all],[100 100]);
imagesc(centers{:}, values.')

colorbar
axis xy 


%%

p=ranksum(inside_bird, outside_bird)

std(inside_bird)
std(outside_bird)

figure
CategoricalScatterplot(([inside_bird;outside_bird]),[ones(length(inside_bird),1);2*ones(length(outside_bird),1)])
 ylim([0 6000])
 
 
 %%
 
 
  %% divide in preferred and not
    
bird_pb=[pb_bird, pb_stim];

categories=[1000*[1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.4,5,6,7,8]]

perc_syll_inside=[];

for ww=1:length(categories)
    [xx,yy]=find(bird_pb(:,2)==categories(ww));
    syllables=sum(yy)
    
    syllables_inside=[];
   
    for www=1:length(xx) 
    if pb_bird(www) < 3652 && pb_bird(www)>1735
        syllables_inside=[ syllables_inside, pb_bird(www)];
    end
end

 
    
    syl_in=length(syllables_inside)
 
    perc=(syl_in/syllables)*100;
    
    perc_syll_inside=[perc_syll_inside, perc];
    
end

[perc_syll_inside(1:4) perc_syll_inside(15:end)]

ranksum(([perc_syll_inside(1:4) perc_syll_inside(15:end)]),(perc_syll_inside(5:14)))
