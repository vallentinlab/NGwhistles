%% Whistle analysis 2.0 for PAIRS matches %%
clear all
close all
clc
% Read all bird's data...
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*.csv');
whistle_table = [];
bird_pitch = [];
stim_pitch = [];
bird_dur = [];
number_bird=[];
number_pb=[];
bad_birds=[2,5,6,7,12,13,14,15,17,20];
good_birds=[1,3,4,8,9,10,11,16,18,19];
bird_id=[];
first_stim=[];
first_resp=[];
%% Load data/ adress folder %%
 for bird=1:20
     
     
if bird == 1
    fprintf("Count with me! :)\n");
end
bird

filename= char(strcat(pathname,'\',filelist(bird,1).name));
[data,text] = xlsread(filename);
%% Read Excel %%
% Only whistle responses (matches)

[idx,idy]=find(data(:,4)==1);

if data(idx(1),1)==0
    id_match=idx(1:2:end);
    id_pb=idx(2:2:end);
elseif data(idx(1),1)==1
    id_match=idx(2:2:end);
    id_pb=idx(1:2:end);
end

bird_w=[];
pb_w=[];

if data(idx(1),1)==0
for ii=1:length(id_pb)-1  
    bi_w=data(id_match(ii):id_pb(ii)-1,5);
    bi_pb=data(id_pb(ii):id_match(ii+1)-1,5);
    mean_pb=nanmean(bi_pb);
    bird_w=[bird_w; bi_w];
    pb_w=[pb_w; mean_pb*ones(length(bi_w),1)];
end

else
    for ii=1:length(id_pb)-1
    bi_pb=data(id_pb(ii):id_match(ii)-1,5);
    bi_w=data(id_match(ii):id_pb(ii+1)-1,5);
 mean_pb=nanmean(bi_pb);
    bird_w=[bird_w; bi_w];
    pb_w=[pb_w; mean_pb*ones(length(bi_w),1)];
end
 end

number=[];
for qq=2:length(data(:,4))-1
    if data(qq,4)>data(qq-1,4) && data(qq,4)>data(qq+1,4)
        number=[number; data(qq,4)];
    end
end
    if data(idx(1),1)==0
number_bird=[number_bird; number(1:2:end-1)];
number_pb=[number_pb; number(2:2:end)];

else
number_pb=[number_pb; number(1:2:end-1)];
number_bird=[number_bird; number(2:2:end)];
    end

bird_pitch = [bird_pitch;bird_w];
stim_pitch = [stim_pitch;pb_w];

for k=1:length(data(:,1))          
if data(k,1)==0  && data(k,4)==1
    bird_freq = data(k,5);
    bird_pitch = [bird_pitch; bird_freq]; % accumulated pitch    
elseif data(k,1)==1 && data(k,4)==1  
    stim_freq = data(k,5); 
    stim_pitch = [stim_pitch; stim_freq];
end
end
% 

[idx,idy]=find(data(:,4)==1);

if data(idx(1),1)==0
    id_match=idx(1:2:end);
    id_pb=idx(2:2:end);
elseif data(idx(1),1)==1
    id_match=idx(2:2:end);
    id_pb=idx(1:2:end);
end

first_stim=[first_stim; (data(id_match,5))];

first_resp=[first_resp; (data(id_pb,5))];

bird_id=[bird_id; bird*ones((length(idx)/2),1)]
% tpbbbb = table(bird_w, pb_w);
% mpbbbb = fitlm(tpbbbb,'bird_w ~ pb_w');
% figure(101)
% ppbbbb = plot(mpbbbb,'color','k',...
%     'Marker','o','MarkerFaceColor', [102 162 162]/256, 'MarkerEdgeColor','none',...
%     'LineWidth',1,'MarkerSize',5);      
%     hold on
%     set(ppbbbb(2),'Color',[195 40 85]/256,'LineWidth',1);
%     ylabel('Focal bird freq (Hz)')
%     xlabel('Stimulus bird freq (Hz)')
%     box off
%     set(gca,'linewidth',1,'FontSize', 14);
%     legend off
%     xlim([0 8000]);
%     ylim([0 8000]);
%     title('');
%     axis square
%     hold on


end  %for loop, all birds

%% Stimuli VS bird frequency (raw) %%
%%



%%
% idd=find(stim_pitch>1932.7948)
% iddd=find(stim_pitch<3425.0649)
% idddd=intersect(idd,iddd)
% 
% stim_pitch=stim_pitch(idddd)
% bird_pitch=bird_pitch(idddd)

tpb = table(first_resp, first_stim);
mpb = fitlm(tpb,'first_resp ~ first_stim')


[rho,pval] = corr(first_resp, first_stim,'Type','Pearson','Rows','complete')


tpb = table(bird_pitch, stim_pitch);
mpb = fitlm(tpb,'bird_pitch ~ stim_pitch')

figure(100)
hold on
ppb = plot(mpb,'color','k',...
    'Marker','o','MarkerFaceColor', [102 162 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(ppb(2),'Color',[195 40 85]/256,'LineWidth',5);
    ylabel('Focal bird freq (Hz)')
    xlabel('Stimulus bird freq (Hz)')
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    xlim([0 8000]);
    ylim([0 8000]);
    title('');
    axis square
    %title('Freq. focal bird VS stimulus bird (raw)');
   
set(gca,'TickDir','out')
          

[rho,pval] = corr(bird_pitch, stim_pitch,'Type','Pearson','Rows','complete')

%% Stimuli VS bird frequency (bins) %%
matrix_data=[stim_pitch, bird_pitch];
sorted_matrix=sortrows(matrix_data,1);
n1=[800:200:8000];

for j=1:length(n1)
    [id_binx,id_biny]=(find(n1(j)<sorted_matrix(:,1) & sorted_matrix(:,1)<n1(j)+200));
    sorted_matrix(id_binx,1)=n1(j);
end



pb_bird=sorted_matrix(:,2);
pb_stim=sorted_matrix(:,1);
tpb = table(pb_bird, pb_stim);
mpb = fitlm(tpb,'pb_bird ~ pb_stim');

figure(2)
ppb =plot(mpb,'color','k','Marker','o','MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
       set(ppb(2), 'Color','k','LineWidth',3);
          ylabel('Focal bird freq (Hz)');
          xlabel('Stimulus bird freq (Hz)');
          box off
          set(gca,'linewidth',1,'FontSize', 14);
          legend off
          xlim([0 8000]);
          ylim([0 8000]);
          title('Freq. focal bird VS stimulus bird (bins)');

%% Density plot %%

data_matrix = nan(10000);
for ww=2:length (n1)
    ids=find(pb_stim==n1(ww));
    data_matrix(1:length(ids),ww)=pb_bird(ids);
 end
       
figure(3)
for qq=1
    subplot(1,length(n1),qq);
    xlim([-10 50]);
    ylim([1000 8000]);
    set(gca,'xtick',[]);
    hold on
end
for qq=2:length(n1);
    subplot(1,length(n1),qq);
    ids=find(pb_stim==n1(qq));
    pb_bird(ids);
    if (sum(pb_stim(ids))>0)==1;
        s1 = scatter_kde(pb_stim(ids)+(10+(rand(length(ids),1))*30),pb_bird(ids));
        s1.SizeData = 15;
        ylim([1000 8000]);    
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(n1(qq)));
        yline(mean(pb_stim(ids)),'-','LineWidth',7); % pb frequency
      % yline(median(pb_bird(ids),'omitnan'),'-','Color', 'r', 'LineWidth',8);
      % yline(mean(pb_bird(ids),'omitnan'),'-','Color', 'm', 'LineWidth',8);

        [gg,bin]=hist(pb_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
      % yline(pb_max,'-','LineWidth',8,'Color', 'r' )
        hold on
        [f, Xi, u] = ksdensity(data_matrix(:,qq));
      % [a,b]=max(f)
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        yline(freq1,'-','LineWidth',7,'Color',[195 40 85]/256) % first peak
% if size(cc,1)>1
%   freq2=Xi(cc(end-1,2));
%   yline(freq2,'-','LineWidth',8,'Color', 'm' ) % second peak
% end
    else 
        ylim([0 8000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(n1(qq)));
        hold on
    end
hold on
end

%% preffered range
q=[0,8000];
qq=[0,8000];

  idd=find(stim_pitch>1735.6);
 iddd=find(stim_pitch<3652.3);

idddd=intersect(idd ,iddd);

 % idddd=[idd ;iddd];

new_stim_pitch=stim_pitch(idddd);
new_bird_pitch=bird_pitch(idddd);

tpb = table(new_bird_pitch, new_stim_pitch);
mpb = fitlm(tpb,'new_bird_pitch ~ new_stim_pitch')

[rho,pval] = corr(new_bird_pitch, new_stim_pitch,'Type','Pearson','Rows','complete')


tpbb = table(bird_pitch, stim_pitch);
mpbb = fitlm(tpbb,'bird_pitch ~ stim_pitch')

figure
hold on
%     plot (q,qq, '--k')
ppb= plot(mpb,'color','k',...
    'Marker','o','MarkerFaceColor', [102 102 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(ppb(2),'Color',[195 140 85]/256,'LineWidth',3);
    ylabel('Focal bird freq (Hz)')
    xlabel('Stimulus bird freq (Hz)')
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    xlim([0 8000]);
    ylim([0 8000]);
    title('');
    xline(1700,'--k');
        xline(3700,'--k');
        set(gca,'DataAspectRatio',[1 1 1])
    %title('Freq. focal bird VS stimulus bird (raw)');
  plot (q,qq, '--k')
set(gca,'TickDir','out')

figure

%  hold on
%     plot (q,qq, '--k')

    xline(1700,'--k');
        xline(3700,'--k');
  
    %title('Freq. focal bird VS stimulus bird (raw)');
    hold on
 ppbb=  plot(mpbb,'color','k',...
    'Marker','o','MarkerFaceColor', [102 162 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(ppbb(2),'Color',[195 40 85]/256,'LineWidth',3);
    ylabel('Focal bird freq (Hz)')
    xlabel('Stimulus bird freq (Hz)')
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    xlim([0 8000]);
    ylim([0 8000]);
    title('');
    axis square
%%

figure
%  hold on
%     plot (q,qq, '--k')

    xline(1700,'--k');
        xline(3700,'--k');
  
    %title('Freq. focal bird VS stimulus bird (raw)');
    hold on
 ppbb=  plot(mpbb,'color','k',...
    'Marker','o','MarkerFaceColor', [102 162 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(ppbb(2),'Color',[195 40 85]/256,'LineWidth',3);
    ylabel('Focal bird freq (Hz)')
    xlabel('Stimulus bird freq (Hz)')
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    xlim([0 8000]);
    ylim([0 8000]);
    title('');
    axis square
    
    hold on
%     plot (q,qq, '--k')
ppb= plot(mpb,'color','k',...
    'Marker','o','MarkerFaceColor', [102 5 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(ppb(2),'Color','r','LineWidth',3);
    ylabel('Focal bird freq (Hz)')
    xlabel('Stimulus bird freq (Hz)')
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    xlim([0 8000]);
    ylim([0 8000]);
    title('');
    xline(1700,'--k');
        xline(3700,'--k');
    %title('Freq. focal bird VS stimulus bird (raw)');
  plot (q,qq, '--k')
  set(gca,'TickDir','out')
    %%
    figure
    raincloud_plot_smooth(bird_pitch, 'line_width', 2,...
    'line_color', [35 70 20]/256, 'bxcl', [35 70 20]/256, 'color', [142 197 69]/256);
   box off
     xlim([0 10000]);
    
    %%
length(bird_pitch)
    
    ttt = table(number_bird, number_pb);
mmm = fitlm(ttt,'number_bird ~ number_pb')

figure
nnn = plot(mmm,'color','k',...
    'Marker','o','MarkerFaceColor', [102 162 162]/256, 'MarkerEdgeColor','none',...
    'LineWidth',1,'MarkerSize',5);      
    hold on
    set(nnn(2),'Color',[195 40 85]/256,'LineWidth',3);
    ylabel('# whistle syllables of matching bird')
    xlabel('# whistle syllables of matched bird')
    box off
    set(gca,'linewidth',1,'FontSize', 14);
    legend off
    title('');
    axis square
    xlim([0 40]);
    ylim([0 40]);
    %title('Freq. focal bird VS stimulus bird (raw)');
   set(gca,'TickDir','out')
    
    
        [R,P] = corrcoef(bird_pitch, stim_pitch,'Rows','complete')  
[rho,pval] = corr(bird_pitch, stim_pitch,'Type','Pearson','Rows','complete')

%%
J=customcolormap([0 0.5 1], [0 0.2 0; 0 0.5 0;1 1 1 ]);

tpb = table(bird_pitch, stim_pitch);
mpb = fitlm(tpb,'bird_pitch ~ stim_pitch')

figure

 [values, centers] = hist3([stim_pitch, bird_pitch ],[75 75]);
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
      xlim([1000 8000]);
    ylim([1000 8000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')
 xline(1735.6,'--k');
        xline(3652.3,'--k');
set(gca,'TickDir','out')
q=[0,9000];
qq=[0,9000];
hold on
  plot (q,qq, ':k')


hold on


%%
figure(8)


tpb = table(bird_pitch,stim_pitch);
mpb = fitlm(tpb,'bird_pitch ~ stim_pitch')
 [values, centers] = hist3([stim_pitch,  bird_pitch ],[75 75]);
imagesc(centers{:}, values.')
colorbar
axis xy 
xlim([1000 8000]);
    ylim([1000 8000]);
        axis square
colormap(J)
hold on
%  ppb = plot(mpb,'color','r','Marker','none');      
%      hold on
%      set(ppb(2), 'Color','r','LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
      xlim([900 9000]);
    ylim([900 9000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')
       %set(gca,'DataAspectRatio',[1 1 1])
    %title('Freq. focal bird VS stimulus bird (raw)');
 xline(1700,'--k');
        xline(3700,'--k');
set(gca,'TickDir','out')
