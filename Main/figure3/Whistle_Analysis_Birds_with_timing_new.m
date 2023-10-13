% Whistle analysis 2.0 %%
clear all
close all
clc
% Read all bird's data...
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*.xlsx');
whistle_table = [];
% Set 1:
pb1_stim = []; pb1_bird = []; pb1_dur = [];
post1_stim = []; post1_bird = []; post1_dur = [];
pre1_stim = []; pre1_bird = []; pre1_dur = []; tempo1=[];wgap1=[];
% Set 2:
pb2_stim = []; pb2_bird = []; pb2_dur = [];
post2_stim = []; post2_bird = []; post2_dur = [];
pre2_stim = []; pre2_bird = []; pre2_dur = []; tempo2=[];wgap2=[];
pb1_stim_end=[];
pb2_stim_end=[];
stim_end1=[];
stim_end2=[];
bird_start1=[];
bird_start2=[];

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

end % for loop, all birds

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
%%
  
%%

figure
subplot(2,1,2)
patch([0 -3 -3 0],[0 0 350 350 ],[0.5 0.5 0.5], 'EdgeColor', 'none')
hold on
for r=1:length(lag)
plot(lag(r),r*1, '.g','MarkerSize', 20)
end 
 xlim([-5 20]);
% xline(-3)
% xline(0)
  box off
  set(gca,'TickDir','out')
  
 %%
 [f, Xi, u] = ksdensity(lag);

 [a]=islocalmin(f);
 [c,b]=findpeaks(f);
 peaks=Xi(b);
 idx=find(a==1);
 throughs=    Xi(idx)   ;
     
 subplot(2,1,1) 
 patch([0 -3 -3 0],[0 0 0.4 0.4 ],[0.5 0.5 0.5], 'EdgeColor', 'none')
hold on
 plot(Xi,f, 'g', 'LineWidth',3)
  xlim([-5 20]);
   ylim([0 0.4]);

 box off
  set(gca,'TickDir','out')
hold on 
 
 %%
  figure
  subplot(2,1,1)
 histogram(lag,15)   
 xlim([-10 20]);
 xline(0)
  box off
  set(gca,'TickDir','out')
 subplot(2,1,2)
 plot(Xi,f)
  xlim([-10 20]);
 xline(0)
 box off
  set(gca,'TickDir','out')
hold on 
plot (throughs,0.05, 'or')
plot (Xi(b),0.15, '*r')

  %%
 idx=find(lag<throughs(1));
 quick_resp_bird=pb_bird(idx);
 quick_resp_stim=pb_stim(idx);
 
 [h,p] = adtest(lag(idx));
 
 
  [ff, Xii, uu] = ksdensity(lag(idx));
  [cc,bb]=findpeaks(ff);
  
  peaks=Xii(bb)
  
 figure
  subplot(2,1,1)
 histogram(lag(idx),'Normalization','probability')   
 xlim([-10 20]);
 xline(0)
  box off
  set(gca,'TickDir','out')
 subplot(2,1,2)
 plot(Xii,ff)
  xlim([-10 20]);
 xline(0)
 box off
  set(gca,'TickDir','out')
hold on 
%%
 
 idx=find(lag>throughs(1));
 
  [h,p] = adtest(lag(idx));

 slow_resp_bird=pb_bird(idx);
 slow_resp_stim=pb_stim(idx);

 
  [ff, Xii, uu] = ksdensity(lag(idx));
  [cc,bb]=findpeaks(ff);
  
  peaks=Xii(bb);
  
 figure
  subplot(2,1,1)
 histogram(lag(idx),'Normalization','probability')   
 xlim([-10 20]);
 xline(0)
  box off
  set(gca,'TickDir','out')
 subplot(2,1,2)
 plot(Xii,ff)
  xlim([-10 20]);
 xline(0)
 box off
  set(gca,'TickDir','out')
hold on 

 %%
% %  pb_bird = pb_bird(~isnan(pb_bird));
% % pb_stim = pb_stim(~isnan(pb_stim));
% 
%  

diff=(quick_resp_bird-quick_resp_stim);

 tpb = table(quick_resp_bird, quick_resp_stim);
mpb = fitlm(tpb,'quick_resp_bird ~ quick_resp_stim')

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
     xlim([0 9000]);
     ylim([0 9000]);
    title ('Freq. bird VS stimuli');
    axis square
 set(gca,'TickDir','out')

 idx=find(lag>throughs(1));
 
 slow_resp_bird=pb_bird(idx);
 slow_resp_stim=pb_stim(idx);
 
diff= slow_resp_bird-slow_resp_stim;
 tpb = table(slow_resp_bird, slow_resp_stim);
mpb = fitlm(tpb,'slow_resp_bird ~ slow_resp_stim')

subplot(1,2,2)
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

 %%
 J=customcolormap([0 0.5 1], [1 1  0 ; 0 0 1;1 1 1 ]);
J=customcolormap([0 0.5 1], [0 0.2 0; 0 0.5 0;1 1 1 ]);

figure
subplot(1,2,1)
 [values, centers] = hist3([ quick_resp_stim quick_resp_bird],[100 100]);
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

  subplot(1,2,2)
 [values, centers] = hist3([ slow_resp_stim slow_resp_bird],[100 100]);
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
