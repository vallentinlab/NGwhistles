clear all; close all; clc

[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*.xlsx');

pb1_stim = []; pb1_bird = []; pb1_dur = []; pb2_stim = []; pb2_bird = []; pb2_dur = [];
stim_end1=[]; stim_end2=[]; bird_start1=[]; bird_start2=[];

all_quick_resp_bird=[];
all_quick_resp_stim=[];
all_slow_resp_bird=[];
all_slow_resp_stim=[];

all_first_peak=[];
all_first_through=[];
all_second_peak=[];
all_lag=[];
all_lag_id=[];
all_lags=nan(500);
all_pb_bird=[ ];
all_pb_stim=[ ];
all_quick_lag=[ ];
all_slow_lag=[ ];

%% Load data/ adress folder %%
for bird=1:11

if bird == 1
    fprintf("Count with me! :)\n");
end
bird

filename= char(strcat(pathname,'\',filelist(bird,1).name));
[~,sheet_name]=xlsfinfo(filename);
data=[];
for k=1:numel(sheet_name)
  [num,txt,raw]=xlsread(filename,sheet_name{k});
%   if isempty(data{1,k})==1
%       data{1,k} = nan;
%   end
data{k}=raw;
dataa{k}=num;
raw=[];
num=[];
end
% Example: x = data{1,2}(:,3);
% data{1,2} will load the 2nd sheet of the 1st file in filelist
% ...(:,3) will get all the rows from the 3rd column of this sheet

%% Read Excel %%
% Damn, this is long. There should be an easier way to do this
% Set 1
% Playback: pb10, pb11, pb12, pb13
if length(data{1,1})==1
    pb10_pitch = nan;    pb10_freq = nan;
    pb10_start = nan;    pb10_end = nan; pb10_stim_end=nan; bird10_start=[]; stim10_end=[];
else
    pb10_pitch =[];   pb10_freq = [];
    bird10_start=[]; stim10_end=[];

     B=convertStringsToChars(data{1,1}(:,4));
        for q=2:length(B)
         if (str2num(B{q}(2:end))==1)==1 | (str2num(B{q}(2:end))==1.1)==1
             
             bird10_start=[bird10_start;(data{1,1}(q,1))];
             stim10_end=[stim10_end; data{1,1}(q,6)];
         
              pb10_pitch =[pb10_pitch; (dataa{1,1}(q,3))]; 
              pb10_freq =[pb10_freq; (dataa{1,1}(q,14))];
         end
        end
    end       

if length(data{1,2})==1
    pb11_pitch = nan;    pb11_freq = nan;
    pb11_start = nan;    pb11_end = nan; pb11_stim_end = nan; bird11_start=[]; stim11_end=[];
else
    pb11_pitch = [];    pb11_freq =[];
    bird11_start=[]; stim11_end=[];

 B=convertStringsToChars(data{1,2}(:,4));
        for q=2:length(B)
         if (str2num(B{q}(2:end))==1)==1 | (str2num(B{q}(2:end))==1.1)==1
             
             bird11_start=[bird11_start;(data{1,2}(q,1))];
             stim11_end=[stim11_end; data{1,2}(q,6)];
             
              pb11_pitch =[pb11_pitch;( dataa{1,2}(q,3))];    pb11_freq =[pb11_freq; (dataa{1,2}(q,14))];
         end
        end
end
if length(data{1,3})==1
    pb12_pitch = nan;    pb12_freq = nan;
    pb12_start = nan;    pb12_end = nan;  pb12_stim_end=nan;  bird12_start=[]; stim12_end=[];
else
    pb12_pitch = [];    pb12_freq = [];
     bird12_start=[]; stim12_end=[];

B=convertStringsToChars(data{1,3}(:,4));
        for q=2:length(B)
         if (str2num(B{q}(2:end))==1)==1 | (str2num(B{q}(2:end))==1.1)==1
             
             bird12_start=[bird12_start;(data{1,3}(q,1))];
             stim12_end=[stim12_end; data{1,3}(q,6)];
             
              pb12_pitch =[pb12_pitch; (dataa{1,3}(q,3))];    pb12_freq =[pb12_freq;(dataa{1,3}(q,14))];
         end
        end
end

if length(data{1,4})==1
    pb13_pitch = nan;    pb13_freq = nan;
    pb13_start = nan;    pb13_end = nan;  pb13_stim_end=nan; bird13_start=[]; stim13_end=[];
else
    pb13_pitch = [];    pb13_freq = [];    
    bird13_start=[]; stim13_end=[];

B=convertStringsToChars(data{1,4}(:,4));
        for q=2:length(B)
         if (str2num(B{q}(2:end))==1)==1 | (str2num(B{q}(2:end))==1.1)==1
             
             bird13_start=[bird13_start;(data{1,4}(q,1))];
             stim13_end=[stim13_end; data{1,4}(q,6)];
             
                 pb13_pitch =[pb13_pitch; (dataa{1,4}(q,3))];    pb13_freq =[pb13_freq; (dataa{1,4}(q,14))];    
         end
        end
end

% Set 2
% Playback: pb20, pb21, pb22
if length(data{1,7})==1
    pb20_pitch = nan;    pb20_freq = nan;
    pb20_start = nan;    pb20_end = nan; pb20_stim_end=nan; bird20_start=[]; stim20_end=[];
else
    pb20_pitch = [];    pb20_freq = [];
    bird20_start=[]; stim20_end=[];
    
B=convertStringsToChars(data{1,7}(:,4));
        for q=2:length(B)
         if (str2num(B{q}(2:end))==1)==1 | (str2num(B{q}(2:end))==1.1)==1
             
             bird20_start=[bird20_start;(data{1,7}(q,1))];
             stim20_end=[stim20_end; data{1,7}(q,6)];
             
                 pb20_pitch =[pb20_pitch; (dataa{1,7}(q,3))];    pb20_freq =[pb20_freq; (dataa{1,7}(q,14))];
         end
        end
end

if length(data{1,8})==1
    pb21_pitch = nan;    pb21_freq = nan;
    pb21_start = nan;    pb21_end = nan; pb21_stim_end=nan;     bird21_start=[]; stim21_end=[];

else
    pb21_pitch = [];    pb21_freq = [];
        bird21_start=[]; stim21_end=[];

B=convertStringsToChars(data{1,8}(:,4));
        for q=2:length(B)
         if (str2num(B{q}(2:end))==1)==1 | (str2num(B{q}(2:end))==1.1)==1
             
             bird21_start=[bird21_start;(data{1,8}(q,1))];
             stim21_end=[stim21_end; data{1,8}(q,6)];
             
                 pb21_pitch =[pb21_pitch; (dataa{1,8}(q,3))];    pb21_freq =[pb21_freq; (dataa{1,8}(q,14))];
         end
        end
end

if length(data{1,9})==1
    pb22_pitch = nan;    pb22_freq = nan;
    pb22_start = nan;    pb22_end = nan; pb22_stim_end=nan;     bird22_start=[]; stim22_end=[];
else
    pb22_pitch = [];    pb22_freq = [];
        bird22_start=[]; stim22_end=[];

B=convertStringsToChars(data{1,9}(:,4));
        for q=2:length(B)
         if (str2num(B{q}(2:end))==1)==1 | (str2num(B{q}(2:end))==1.1)==1
             
             bird22_start=[bird22_start;(data{1,9}(q,1))];
             stim22_end=[stim22_end; data{1,9}(q,6)];
             
                 pb22_pitch =[pb22_pitch; (dataa{1,9}(q,3))];    pb22_freq =[pb22_freq; (dataa{1,9}(q,14))];
         end
        end
end
% Storing data from all birds in new variables
pb1_stim = [pb1_stim; pb10_freq; pb11_freq; pb12_freq; pb13_freq];
pb1_bird = [pb1_bird; pb10_pitch; pb11_pitch; pb12_pitch; pb13_pitch];
pb2_stim = [pb2_stim; pb20_freq; pb21_freq; pb22_freq];
pb2_bird = [pb2_bird; pb20_pitch; pb21_pitch; pb22_pitch];
stim_end1=[stim_end1;stim10_end;stim11_end;stim12_end;stim13_end];
stim_end2=[stim_end2;stim20_end;stim21_end;stim22_end];
bird_start1=[bird_start1;bird10_start;bird11_start;bird12_start;bird13_start];
bird_start2=[bird_start2;bird20_start;bird21_start;bird22_start];

pb_bird=[pb1_bird; pb2_bird];
pb_stim=[pb1_stim; pb2_stim];
bird_start=[bird_start1;bird_start2];
stim_end=[stim_end1;stim_end2];
bird_start=cell2mat(bird_start);
stim_end=cell2mat(stim_end);
lag=bird_start-stim_end;
lag=lag/44100;
bird_start=bird_start/44100;
stim_end=stim_end/44100;

[f, Xi, u] = ksdensity(lag);
 [a]=islocalmin(f);
 idx=find(a==1);
 throughs=    Xi(idx);
 idx=find(lag<throughs(1));
 quick_resp_bird=pb_bird(idx);
 quick_resp_stim=pb_stim(idx);
 
 
 idx=find(lag>throughs(1));
 slow_resp_bird=pb_bird(idx);
 slow_resp_stim=pb_stim(idx);
 
e=[rand rand rand];
 
diff=(quick_resp_bird-quick_resp_stim);
 
tpb = table(quick_resp_bird, quick_resp_stim);
mpb = fitlm(tpb,'quick_resp_bird ~ quick_resp_stim')

[rho,pval] = corr(quick_resp_bird, quick_resp_stim,'Type','Pearson','Rows','complete')

figure(1)
subplot(2,2,3)
if pval<0.05
ppb = plot(mpb,'color',e,'Marker','none','MarkerSize',4,'MarkerFaceColor', e, 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color',e,'LineWidth',3);
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

end
hold on
%  
diff=(slow_resp_bird-slow_resp_stim);
 tpb = table(slow_resp_bird, slow_resp_stim);
mpb = fitlm(tpb,'slow_resp_bird ~ slow_resp_stim')

[rho,pval] = corr(slow_resp_bird, slow_resp_stim,'Type','Pearson','Rows','complete')


subplot(2,2,4)
if pval<0.05
ppb = plot(mpb,'color',e,'Marker','none','MarkerSize',4,'MarkerFaceColor', e, 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color',e,'LineWidth',3);
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
end
 hold on
 
 
  [f, Xi, u] = ksdensity(lag);

 
subplot(2,2,1:2)
 %patch([0 -3 -3 0],[0 0 0.4 0.4 ],[0.5 0.5 0.5], 'EdgeColor', 'none')
hold on
 plot(Xi,f, 'color', e, 'LineWidth',2)
  xlim([-5 20]);
   ylim([0 0.5]);

 box off
  set(gca,'TickDir','out')
hold on 
 

 figure
subplot(2,1,2)
patch([0 -3 -3 0],[0 0 350 350 ],[0.5 0.5 0.5], 'EdgeColor', 'none')
hold on
for r=1:length(lag)
plot(lag(r),r*1, '.', 'color', e,'MarkerSize', 20)
end 
 xlim([-5 20]);
% xline(-3)
% xline(0)
  box off
  set(gca,'TickDir','out')

 subplot(2,1,1) 
 patch([0 -3 -3 0],[0 0 0.4 0.4 ],[0.5 0.5 0.5], 'EdgeColor', 'none')
hold on
 plot(Xi,f, 'color', e, 'LineWidth',3)
  xlim([-5 20]);
   ylim([0 0.4]);

 box off
  set(gca,'TickDir','out')
hold on 

 
 
[aa,bb]=findpeaks(f);
       
        peak1=Xi(bb(1));
        peak2=Xi(bb(2));


 
all_first_peak=[all_first_peak;    peak1       ];
all_first_through=[all_first_through;    throughs(1)   ];
all_second_peak=[all_second_peak;   peak2       ];

all_quick_resp_bird=[all_quick_resp_bird;quick_resp_bird];
all_quick_resp_stim=[all_quick_resp_stim;quick_resp_stim];

quick_lag=lag(lag<throughs(1))

all_quick_lag=[all_quick_lag; quick_lag];
 
all_slow_resp_bird=[all_slow_resp_bird;slow_resp_bird];
all_slow_resp_stim=[all_slow_resp_stim;slow_resp_stim];

slow_lag=lag(lag>throughs(1))

all_slow_lag=[all_slow_lag; slow_lag];

all_pb_bird=[all_pb_bird; pb_bird];
all_pb_stim=[all_pb_stim; pb_stim];
all_lags(1:length(lag)',bird)=lag';

all_lag=[all_lag;lag];
all_lag_id=[all_lag_id;bird*ones(length(lag),1)];

pb1_stim = [];
pb1_bird = [];
pb2_stim = [];
pb2_bird = [];
stim_end1=[];
stim_end2=[];
bird_start1=[];
bird_start2=[];

pb_bird=[];
pb_stim=[];
bird_start=[];
stim_end=[];
bird_start=[];
stim_end=[];
lag=[];
end % for loop, all birds
%%
[ooo,aaa]=find(all_quick_lag<0);
ov_lag=all_quick_lag(ooo);

[ooo,aaa]=find(all_quick_lag>0);
alt_lag=all_quick_lag(ooo);

jj=(length(all_quick_lag)/length(all_lag))*100;

jjj=(length(all_slow_lag)/length(all_lag))*100;

figure
subplot(1,2,1)
pie([jj jjj])
subplot(1,2,2)
j=(length(alt_lag)/length(all_quick_lag))*100;

jj=(length(ov_lag)/length(all_quick_lag))*100;
pie([j jj])

%%

n = all_slow_resp_stim(~isnan(all_slow_resp_stim));


[uuu,eee]=find(n<1700);
[uuuu,eeee]=find(n>3700);
 uu=[uuu; uuuu];%    intersect(uuu, uuuu);  %

outside_quick=n(uu);

hh=(length(outside_quick)/length(n))*100



%%
[r,t]=find(all_quick_resp_stim<1600);
[rr,tt]=find(all_quick_resp_stim>3800);


[ir,it]=find(all_quick_resp_stim>1600);
[irr,itt]=find(all_quick_resp_stim<3800);

rrr=[r;rr];

irrr=intersect(ir,irr);

central_quick_bird=all_quick_resp_bird(irrr);
central_quick_stim=all_quick_resp_stim(irrr);
central_quick_lag=all_quick_lag(irrr);

central_quick_diff=abs(central_quick_bird-central_quick_stim);


outcentral_quick_bird=all_quick_resp_bird(rrr);
outcentral_quick_stim=all_quick_resp_stim(rrr);
outcentral_quick_lag=all_quick_lag(rrr);

outcentral_quick_diff=abs(outcentral_quick_bird-outcentral_quick_stim);

[rq,tq]=find(all_slow_resp_stim<1600);
[rrq,ttq]=find(all_slow_resp_stim>3800);
rrrq=intersect(rq,rrq);

rrrq=[rq;rrq];



outcentral_slow_bird=all_slow_resp_bird(rrrq);
outcentral_slow_stim=all_slow_resp_stim(rrrq);
outcentral_slow_lag=all_slow_lag(rrrq);

outcentral_slow_diff=abs(outcentral_slow_bird-outcentral_slow_stim);

outcentral_stim=[outcentral_quick_stim;outcentral_slow_stim]
outcentral_diff=[outcentral_quick_diff;outcentral_slow_diff];
outcentral_lag=[outcentral_quick_lag;outcentral_slow_lag];


 tpb = table(outcentral_quick_lag,outcentral_quick_stim);
mpb = fitlm(tpb,'outcentral_quick_lag~ outcentral_quick_stim ')

figure
ppb = plot(mpb,'color','r','Marker','o','MarkerSize',4,'MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','r','LineWidth',3);



%%
 
 J=customcolormap([0 0.5 1], [0 0.2 0; 0 0.5 0;1 1 1 ]);

figure
hold on
 [values, centers] = hist3([ all_lag, all_lag_id],[100 21]);
imagesc(centers{:}, values.');
colorbar
axis xy 
xlim([-5 20]);
    ylim([0 12]);
        axis square
colormap(J)
xline(0, '--')
xline(-3, '--')
box off
set(gca,'TickDir','out')





%%
figure  
patch([0 -3 -3 0],[0 0 7000 7000 ],[0.5 0.5 0.5], 'EdgeColor', 'none')
 violinplot(([all_first_peak;all_first_through;all_second_peak]),...
     [ones(length(all_first_peak),1);2*ones(length(all_first_through),1);3*ones(length(all_second_peak),1)])
box off
 ylim([-5 10]);
axis square
 set(gca,'TickDir','out')

%%

[f,ff]=find(all_quick_lag<0)

diff1=abs(all_quick_resp_bird(f)-all_quick_resp_stim(f));

fff=all_quick_lag(f);

 tpb = table(fff,diff1);
mpb = fitlm(tpb,'diff1~ fff ')

figure
subplot(1,2,1)
ppb = plot(mpb,'color','r','Marker','o','MarkerSize',4,'MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','r','LineWidth',3);
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

 idx=find(lag>throughs(1));
 
 slow_resp_bird=pb_bird(idx);
 slow_resp_stim=pb_stim(idx);
 
%  
diff2=(all_slow_resp_bird-all_slow_resp_stim);

 tpb = table(diff2, all_slow_lag);
mpb = fitlm(tpb,'diff2 ~ all_slow_lag')

subplot(1,2,2)
ppb = plot(mpb,'color','r','Marker','o','MarkerSize',4,'MarkerFaceColor', 'r', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','r','LineWidth',3);
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
 
%%
% l=[all_quick_lag;all_slow_lag];
% d=(abs([diff1;diff2]));
% 
% 
%  J=customcolormap([0 0.5 1], [1 1  0 ; 0 0 1;1 1 1 ]);
% J=customcolormap([0 0.5 1], [0 0.2 0; 0 0.5 0;1 1 1 ]);

figure
subplot(1,2,1)
%  [values, centers] = hist3([l d  ],[100 100]);
% imagesc(centers{:}, values.')
% 
% colorbar
% colormap(J)
% axis xy 
% xlim([0 8000]);
%     ylim([0 8000]);
%         axis square
% 
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
  yline(1735,'--')
 yline(3652,'--')

q=[0,9000];
qq=[0,9000];
hold on
  plot (q,qq, ':k')

  subplot(1,2,2)
 [values, centers] = hist3([ all_slow_resp_stim all_slow_resp_bird],[100 100]);
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
  yline(1735,'--')
 yline(3652,'--')

q=[0,9000];
qq=[0,9000];
hold on
  plot (q,qq, ':k')
%%

diff1=((all_quick_resp_bird)-(all_quick_resp_stim));
diff2=((all_slow_resp_bird)-(all_slow_resp_stim));


d=(abs([diff1;diff2]));
l=[all_quick_lag;all_slow_lag];
 tpb = table(d, l);
mpb = fitlm(tpb,'d ~ l');

figure
scatter(all_quick_lag,abs(diff1));
hold on
scatter(all_slow_lag,abs(diff2));

[rho,pval] = corr( l,d,'Type','Pearson','Rows','complete')


figure
ppb = plot(mpb,'color','k','Marker','none','MarkerSize',4,'MarkerFaceColor', 'none', 'MarkerEdgeColor','none');      
    hold on
    set(ppb(2), 'Color','k','LineWidth',3);
    ylabel('Bird freq (Hz)');
    xlabel('Stimulus freq (Hz)');
    box off
    set(gca,'linewidth',1,'FontSize', 14)
    legend off
     xlim([-5 20]);
     ylim([0 7000]);
    title ('Freq. bird VS stimuli');
    set(gca,'TickDir','out')


group1=1*ones(length(diff1),1);
group2=2*ones(length(diff2),1);
groupp=[group1; group2];

figure
beeswarm(groupp,abs([diff1;diff2]))
 xlim([-5 10]);
figure
boxplot(abs([diff1;diff2]),groupp)
figure
CategoricalScatterplot(abs([diff1;diff2]),groupp)
anova1(abs([diff1;diff2]), groupp)

[p,h]=ranksum(abs(diff1),abs(diff2))

3-(-all_first_peak)
%%

[o,oo]=find(all_quick_resp_stim>1700 & all_quick_resp_stim<3700 );

(sum(oo)/length(all_quick_resp_stim))*100

[e,ee]=find(all_slow_resp_stim>1700 & all_slow_resp_stim<3700 );

CC=(sum(ee)/length(all_slow_resp_stim))*100
C=(sum(oo)/length(all_quick_resp_stim))*100

labels = {'to PBs central','to PBs outside'};
figure
subplot(1,2,1)
pie([C 100-C], labels)
title('Early responses')
subplot(1,2,2)
pie([CC 100-CC], labels)
title('Late responses')

aaa=table([48; 52], [46 ;54])
[h,p,stats] = fishertest(aaa)

