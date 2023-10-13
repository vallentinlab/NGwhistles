clear all;close all;clc


% Read all bird's data...
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*.csv');
whistle_table = [];

%% Load data/ adress folder %%

  all_playback=[];
    all_response=[];
X=[]; Y=[];
playback=[];
response=[];
w_id=[];
for bird=1:10

filename= char(strcat(pathname,'\',filelist(bird,1).name));
[~,sheet_name]=xlsfinfo(filename);
  data=xlsread(filename);
  
  W=data(:,1);
  
 [ firs_syl, q]=find(data(:,1)==1);


playback=[playback; data(firs_syl,12)];
  response=[response; data(firs_syl,6)];
  w_id=[w_id; data(:,1)];
  
  all_playback=[all_playback; data(:,12)];
      
  all_response=[all_response; data(:,6)];
  
  ppp=[2400 2800 3200 5000 6000];
  e= [rand rand rand];
  figure(1)
  scatter(ppp,ppp,"_", 'MarkerEdgeColor', 'k' ,'MarkerFaceColor', 'k','LineWidth',5)
hold on
  scatter(data(firs_syl,12)+(10+(rand(length(data(firs_syl,12)),1))*75),data(firs_syl,6),50,'MarkerEdgeColor', 'none' ,'MarkerFaceColor',e )
  hold on
       xlim([1000 7000]);
   ylim([1000 7000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')
  
  
  
 figure(2)
 
 sti=data(firs_syl,12);
 resp=data(firs_syl,6);
 
 tpb = table(sti, resp);
mpb = fitlm(tpb,'sti ~ resp')

     legend off
     xlim([1000 7000]);
   ylim([1000 7000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')
     hold on
ppb = plot(mpb,'color',e,'Marker','o','MarkerFaceColor', 'none', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color',e,'LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off

 
 

end

%%
% playback=playback(x);
% response=response(x);

pb_stim=[response playback]

figure
scatter(playback,response)


delta=playback-response

r=100*[24 28 32 50 60];

pp=[];

for t=1:5
    [o,oo]=find(all_playback==r(t));
   
    [u,uu]=find(all_response(o)< r(t)+100 & all_response(o)>r(t)-100);
    
    pp=[pp, (sum(uu)/sum(oo))*100];
 
end



%%
PB_shuff=all_playback(randperm(length(all_playback)));

RESP_shuff=all_response(randperm(length(all_response)));

shuf_pp=[];

for t=1:5
    [o,oo]=find(PB_shuff==r(t));
   
    [u,uu]=find(RESP_shuff(o)< r(t)+100 & RESP_shuff(o)>r(t)-100);
    
    shuf_pp=[shuf_pp, (sum(uu)/sum(oo))*100];
 
end


%%

figure
plot([2*ones(1,length(pp))], (pp),...
    'o')
hold on
plot([1*ones(1,length(shuf_pp))], (shuf_pp),...
    'o')
hold on
 plot([1*ones(1,length(shuf_pp));2* ones(1,length(pp))],[shuf_pp;pp],...
     '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on

box off
    xlim([0 3]);
 
        axis square
 set(gca,'TickDir','out')
 
[p,h]=ttest(pp,shuf_pp)

mean(pp)

%%
    PB_shuff=playback(randperm(length(playback)));

RESP_shuff=response(randperm(length(response)));

shuff_pb_stim=[RESP_shuff PB_shuff]

figure
scatter(PB_shuff,RESP_shuff)

Shuf_delta=PB_shuff-RESP_shuff;


[rho,pval] = corr(RESP_shuff, PB_shuff,'Type','Pearson','Rows','complete')

%%
figure
histogram(delta, 100)
xline(median(delta))


    
     %%
figure
raincloud_plot((delta))
[p,h]=ttest(delta)
  xlim([-5000 5000]);
box off
    set(gca,'linewidth',1,'FontSize', 14)
        axis square
 set(gca,'TickDir','out')
 

 %%
 jj = [1000:200:7000];

data_matrix = nan(200);

 for ww=1:length (jj)
  ids=find(playback==jj(ww));
  data_matrix(1:length(ids),ww)=response(ids);
 end
 
 
h=[]
p=[];
 for ww=1:length (jj)
 if (sum(data_matrix(:,ww),'omitnan')>0)==1
     [pp,hh]=signtest((data_matrix(:,ww)),jj(ww) )
     p=[p,pp]
     h=[h,hh]
 else
     p=[p,0]
     h=[h,0]
 end
 end
 
 
 
figure
plot([1, 31],[1000, 7000], '--')
hold on
violinplot(data_matrix);
box off
set(gca,'TickDir','out')
yline(0, '-')
axis square
ylim([1000 7000])
xlim([1 31])
xticks([8 10 12 21 26])
xticklabels({'2400','2800','3200','5000','6000'})


%%

medians=[];
for q=1:length (jj)
mmm=median(data_matrix(:,q),'omitnan')
medians=[medians, mmm];
end

medians=medians(~isnan(medians));
stim=[2400, 2800, 3200, 5000, 6000];

[p,h]=signrank(medians, stim)

diff=medians-stim

[p,h]=signrank(diff)

figure
CategoricalScatterplot(diff',ones(length(diff),1))
xlim([0.5 1.5])
box off
set(gca,'TickDir','out')
yline(0, '-')
axis square

%%

figure
plot([2*ones(1,length(medians))], (medians),...
    'o','MarkerFaceColor', 'r', 'MarkerEdgeColor','r', 'MarkerSize', 7)
hold on
plot([1*ones(1,length(stim))], (stim),...
    'o','MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 7)
hold on
plot([1*ones(1,length(stim));2* ones(1,length(stim))],[stim; medians],...
    '-', 'Color', [0.5 0.5 0.5],'LineWidth',1)
hold on
xlim([0 3])
ylim([1000 7000])
box off
xticks([1 2])
xticklabels({'Playback','Response'})
yticks([ 2000 3000 4000 5000 6000])
%legend(sprintf('p = %0.3f',p_anti))
%title('Percentage of whistle songs')
axis square
hold on
set(gca,'TickDir','out')

%% shuffle

mmm=median(data_matrix,'omitnan');
medians=mmm(~isnan(mmm))


M=[];

for r=1:1000
shuf_data_matrix = nan(200);

PB_shuff=playback(randperm(length(playback)));

RESP_shuff=response(randperm(length(response)));


 for ww=1:length (jj)
  ids=find(PB_shuff==jj(ww));
  shuf_data_matrix(1:length(ids),ww)=RESP_shuff(ids);
 end
 
m=median(shuf_data_matrix,'omitnan');
shuf_medians=m(~isnan(m));

M=[M; shuf_medians];

end


%%

figure

for h=1:length(stim)
%       scatter(stim,stim, 200, '_g',"LineWidth",5)
    scatter(stim(h)*ones(height(M),1),M(:,h),"_",'k')
    hold on
      
  [e,ee]=find(M(:,h)>medians(h));
  sum(ee)
  
    scatter(stim,medians, 100, '*r')
  
end

axis square
    box off
     xlim([1000 7000]);
   ylim([1000 7000]);
 set(gca,'TickDir','out')
  
 xticks([2400 2800 3200 5000 6000 ])
plot([1, 7000],[1, 7000], 'k:')

%%


figure

subplot(1,2,1)

for qww=1:length(medians)
   e= [rand rand rand];
stem(qww,medians(qww)-stim(qww),'filled','LineWidth',2,'LineStyle','-.','Color',[e],...
     'MarkerFaceColor',[e], 'MarkerSize', 8)
 hold on
end
xlim([0 6])
 ylim([-1500 1500])
 
%axis square

box off
set(gca,'TickDir','out')
% set(gca, 'yScale', 'log')

subplot(1,2,2)
CategoricalScatterplot(diff',ones(length(diff),1))
xlim([0.5 1.5])
box off
set(gca,'TickDir','out')
yline(0, '-')
axis square
 ylim([-1500 1500])

 mean(diff)
 std(diff)
[p,h]= ttest(diff)
 
 %%
 

real_matrix=[80.6 16.7 0 0 2.8; 27.3 45.5 27.3 0 0; 0 22.2  77.8 0 0; 0 11.1 33.3 22.2 33.3; 0 0 7.7 0 92.3];
shuffled_matrix=[72.2 2.8 25.0 0 0; 59.1 4.5 36.4 0 0; 59.1 11.1 37.0 0 0; 55.6 0 44.4 0 0; 46.2 0 53.8 0 0]

%%

figure
imagesc(real_matrix)
cmap=(summer(2048));
colormap(cmap);
cmap(1,:)=[1,1,1];
colormap(cmap)
colorbar
box off
box off
axis square
set(gca,'TickDir','out')
set(gca,'color','none')
caxis([0 100])

figure
imagesc(shuffled_matrix)
cmap=(summer(2048));
colormap(cmap);
cmap(1,:)=[1,1,1];
colormap(cmap)
colorbar
box off
box off
axis square
set(gca,'TickDir','out')
set(gca,'color','none')
caxis([0 100])


%%

mean_diag=[];

ww=[-4:1:4];

for qq=1:9

x2 = diag(real_matrix,ww(qq));
mean_diag=[mean_diag, mean(x1)];

end


figure
subplot(1,2,1)
plot(ww,mean_diag)
box off
box off
axis square
set(gca,'TickDir','out')
xlim([-4 4])
ylim([0 75])



mean_diag=[];

ww=[-4:1:4];

for qq=1:9

x1 = diag(shuffled_matrix,ww(qq));
mean_diag=[mean_diag, mean(x1)];

end


subplot(1,2,2)
plot(ww,mean_diag)
box off
box off
axis square
set(gca,'TickDir','out')

xlim([-4 4])
ylim([0 75])



%%

%%
tpb = table(playback, response);
mpb = fitlm(tpb,'response ~ playback')

[rho,pval] = corr(response, playback,'Type','Pearson','Rows','complete')

    
J=customcolormap([0 0.5 1], [0 0.2 0; 0 0.5 0;1 1 1 ]);

figure
 [values, centers] = hist3([   playback response],[60 60]);
imagesc(centers{:}, values.')

colorbar
axis xy 
axis square
colormap(J)
hold on
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
     legend off
     xlim([1000 7000]);
   ylim([1000 7000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')
     hold on
ppb = plot(mpb,'color','b','Marker','o','MarkerFaceColor', 'none', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','k','LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off

    xticks([2400, 2800, 3200 5000, 6000])

 hold on
 q=[0,10000];
 plot(q,q, '--')
    plot([2400 2800 3200 5000 6000],medians, 'r_')
%%




ww=[2400 2800 3200 5000 6000];

data_mat=nan(10000);

for w=1:length (ww)
    ids=find(all_playback==ww(w));
    data_mat(1:length(ids),w)=all_response(ids);
end


figure

STDs=[];

peak_inside=[];
accurate=[];

for qq=1:length(ww);
    subplot(1,length(ww),qq);
    ids=find(all_playback==ww(qq));
    all_response(ids);
    zzz=all_response(ids);
    ww(qq);
    %[p,h]=ranksum(zzz,pre_bird);
        scatter(all_playback(ids)+(10+(rand(length(ids),1))*30),all_response(ids));
        ylim([1000 7000]);
        set(gca,'xtick',[]);
        xlabel(convertCharsToStrings(ww(qq)));
      hold on
       set(gca,'TickDir','out')
       
       yline(ww(qq), 'k_')
              STDs=[STDs ; std(all_response(ids))]

             [a,aa]= find(all_response(ids)<ww(qq)+100 & all_response(ids)>ww(qq)-100);
              
             accurate=[accurate, (sum(aa)/length(all_response(ids)))*100];
              
end
%%

accurate_shuff=[];
figure
for qq=1:length(ww);
    subplot(1,length(ww),qq);
    ids=find(PB_shuff==ww(qq));
    RESP_shuff(ids);
    zzz=RESP_shuff(ids);
    ww(qq);
    %[p,h]=ranksum(zzz,pre_bird);
        scatter(PB_shuff(ids)+(10+(rand(length(ids),1))*30),RESP_shuff(ids));
        ylim([1000 7000]);
        set(gca,'xtick',[]);
        xlabel(convertCharsToStrings(ww(qq)));
      hold on
       set(gca,'TickDir','out')
       
       yline(ww(qq), 'k_')
              STDs=[STDs ; std(RESP_shuff(ids))]

             [a,aa]= find(RESP_shuff(ids)<ww(qq)+100 & RESP_shuff(ids)>ww(qq)-100);
              
             accurate_shuff=[accurate_shuff, (sum(aa)/length(RESP_shuff(ids)))*100];
              
end












%%
 tpb = table(playback, response);
mpb = fitlm(tpb,'playback ~ response')

[rho,pval] = corr(response, playback,'Type','Pearson','Rows','complete')


%%
%%
x=[31 33 2 14 7 3 2 10 2 3];
y=[10 15 2 7 3 0 2 0 0  0];

w=x+y;

per=(x./w)*100;
figure
scatter((1+0.3*(rand(length(per),1))),per )
xlim([0 2])
ylim([0 100])
mean(per)