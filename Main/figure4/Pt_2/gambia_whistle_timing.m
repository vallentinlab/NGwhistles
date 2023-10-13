clear all;close all;clc


% Read all bird's data...
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*.csv');
whistle_table = [];

%% Load data/ adress folder %%

  all_playback=[];
    all_response=[];

playback=[];
response=[];
off_pb=[];
on_bird=[];

for bird=1:10

filename= char(strcat(pathname,'\',filelist(bird,1).name));
[~,sheet_name]=xlsfinfo(filename);
  data=xlsread(filename);
  
 
  
  playback=[playback; data(:,3)];
  response=[response; data(:,6)];
 
  off_pb=[off_pb; data(:,2)];
  on_bird=[on_bird; data(:,4)];
  
  lag_OFF=data(:,4)-data(:,2);
  
  figure(1)
 subplot(1,4,1)
  plot(lag_OFF,bird*1, '|', 'color', [0 0.5 0],'MarkerSize',15, 'LineWidth', 3)
hold on
 box off
  set(gca,'TickDir','out')
   xlim([-5 15]);
      ylim([0 11]);

   axis square
xline(0, '--')
box off
set(gca,'TickDir','out')

ylabel('Bird #')

xlabel('latency to whistle song offset (sec)')
  hold on
  
end

%%

error=response-playback;
latency=on_bird-off_pb;
       
      figure(1)
 subplot(1,4,2)
     histogram(latency,30, 'Normalization', 'probability')

hold on
  
     [f, Xi, u] = ksdensity(latency);
    hold on
    xline(0, '--')
  
       plot(Xi, f)
 [a,b]=findpeaks(f);

        c=[a;b];
        cc=sortrows(c');
       
   xlim([-5 15]);
   ylim([0 0.42]);
   
   box off
       axis square
 set(gca,'TickDir','out')
   ylabel('Propability')

xlabel('latency to whistle song offset (sec)')


 [c]=islocalmin(f);
 idx=find(c==1);
 through=    Xi(idx);
 
  
  %%
  
  err=abs(error);
 
  [i,ii]=find(latency<through(1));

    [a,aa]=find(latency>through(1));

  ranksum(err(i), err(a))
    
median(err(i))
 
  figure(1)
 subplot(1,4,3)
scatter( latency(i),err(i))
hold on
scatter( latency(a),err(a))
 ylim([0 5000]);

  legend off
   xlabel('latency to whistle song offset (sec)')
    xline(0, '--')
   ylabel('pitch error (Hz)')
   xlim([-5 15]);
  axis square
 set(gca,'TickDir','out')
  
group1=1*ones(length(i),1);
group2=2*ones(length(a),1);
groupp=[group1; group2];


  figure(1)
 subplot(1,4,4)
CategoricalScatterplot(abs([err(i);err(a)]),groupp)
axis square
 set(gca,'TickDir','out')
 ylim([0 5000]);

%%
[g,gg]=find(latency<through(1));
[ge,gge]=find(latency>through(1));


figure (11)
subplot(1,2,1)

pie([sum(gg) sum(gge)])


[g,gg]=find(latency<through(1) & latency<0);
[ge,gge]=find(latency<through(1) & latency>0);


figure(11)
subplot(1,2,2)
pie([sum(gg) sum(gge)])



