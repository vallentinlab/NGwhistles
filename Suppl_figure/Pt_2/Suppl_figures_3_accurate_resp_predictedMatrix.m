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

bird_id1_stim=[];
bird_id2_stim=[];

bird_id1_pre=[];
bird_id2_pre=[];

bird_id1_post=[];
bird_id2_post=[];


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
    pb20_pitch = [];    pb20_freq = [];
    pb20_start = [];    pb20_end = [];
else
    pb20_pitch = data{1,7}(:,3);    pb20_freq = data{1,7}(:,14);
    pb20_start = data{1,7}(:,1);    pb20_end = data{1,7}(:,2);
end
if isempty(data{1,8})==1
    pb21_pitch = [];    pb21_freq = [];
    pb21_start = [];    pb21_end = [];
else
    pb21_pitch = data{1,8}(:,3);    pb21_freq = data{1,8}(:,14);
    pb21_start = data{1,8}(:,1);    pb21_end = data{1,8}(:,2);
end
if isempty(data{1,9})==1
    pb22_pitch = [];    pb22_freq = [];
    pb22_start = [];    pb22_end = [];
else
    pb22_pitch = data{1,9}(:,3);    pb22_freq = data{1,9}(:,14);
    pb22_start = data{1,9}(:,1);    pb22_end = data{1,9}(:,2);
end

% Postcontrol: post2
if isempty(data{1,10})==1
    post2_pitch = [];    post2_freq = [];
    post2_start = [];    post2_end = [];
else
    post2_pitch = data{1,10}(:,3);    post2_freq = data{1,10}(:,14);
    post2_start = data{1,10}(:,1);    post2_end = data{1,10}(:,2);
end

% Precontrol: pre20, pre21
if isempty(data{1,11})==1
    pre20_pitch = [];    pre20_freq = [];
    pre20_start = [];    pre20_end = [];
else
    pre20_pitch = data{1,11}(:,3);    pre20_freq = data{1,11}(:,14);
    pre20_start = data{1,11}(:,1);    pre20_end = data{1,11}(:,2);
end
if isempty(data{1,12})==1
    pre21_pitch = [];    pre21_freq = [];
    pre21_start = [];    pre21_end = [];
else
    pre21_pitch = data{1,12}(:,3);    pre21_freq = data{1,12}(:,14);
    pre21_start = data{1,12}(:,1);    pre21_end = data{1,12}(:,2);
end

% Storing data from all birds in new variables
pb1_stim = [pb1_stim; pb10_freq; pb11_freq; pb12_freq; pb13_freq];
pb1_bird = [pb1_bird; pb10_pitch; pb11_pitch; pb12_pitch; pb13_pitch];
pb2_stim = [pb2_stim; pb20_freq; pb21_freq; pb22_freq];
pb2_bird = [pb2_bird; pb20_pitch; pb21_pitch; pb22_pitch];

post1_stim = [post1_stim; post1_freq];
post1_bird = [post1_bird; post1_pitch];
post2_stim = [post2_stim; post2_freq];
post2_bird = [post2_bird; post2_pitch];

pre1_stim = [pre1_stim; pre1_freq];
pre1_bird = [pre1_bird; pre1_pitch];
pre2_stim = [pre2_stim; pre20_freq; pre21_freq];
pre2_bird = [pre2_bird; pre20_pitch; pre21_pitch];


bird_id1_stim=[bird_id1_stim zeros*(1:length(pb10_freq))+bird zeros*(1:length(pb11_freq))+bird ...
    zeros*(1:length(pb12_freq))+bird  zeros*(1:length(pb13_freq))+bird];
bird_id2_stim=[bird_id2_stim zeros*(1:length(pb20_freq))+bird zeros*(1:length(pb21_freq))+bird zeros*(1:length(pb22_freq))+bird];

bird_id1_pre=[bird_id1_pre zeros*(1:length(pre1_pitch))+bird];
bird_id2_pre=[bird_id2_pre zeros*(1:length(pre20_pitch))+bird  zeros*(1:length(pre21_pitch))+bird];

bird_id1_post=[bird_id1_post zeros*(1:length(post1_pitch))+bird];
bird_id2_post=[bird_id2_post zeros*(1:length(post2_pitch))+bird];



end % for loop, all birds


id_pre=[bird_id1_pre bird_id2_pre]';
id_post=[bird_id1_post bird_id2_post]';
id_pb=[bird_id1_stim bird_id2_stim]';

post_bird=[post1_bird; post2_bird];
pre_bird=[pre1_bird; pre2_bird];
pb_stim=[pb1_stim; pb2_stim];
pb_bird=[pb1_bird;pb2_bird];

 clearvars -except post_bird pre_bird pb_stim pb_bird post1_bird post2_bird ... 
     pre1_bird pre2_bird pb1_stim pb2_stim pb1_bird pb2_bird id_pre id_post id_pb

bird_stim_id=[pb_bird, pb_stim];

[u,uu]=find(pb_stim>1700 & pb_stim<3800);


stim_resp=[pb_bird(u), pb_stim(u)];

shuff_bird=pb_bird(randperm(length(pb_bird)));
shuff_stim=pb_stim(randperm(length(pb_stim)));

[u,uu]=find(shuff_stim>1700 & shuff_stim<3800);



shuff=[shuff_bird, shuff_stim];

pre=[pre_bird, id_pre];
%% calculate Perc % outside for PBs inside and outside

limit_low=1735;
limit_high=3652;

jj = [600:200:8600];

data_matrix = nan(10000);

data_matrix(1:length(pre_bird),1) = pre_bird;
data_matrix(1:length(post_bird),41) = post_bird;

 for ww=2:length (jj)
     
  ids=find(bird_stim_id(:,2)==jj(ww));
  
  x=bird_stim_id(:,1);
  
  data_matrix(1:length(ids),ww)=x(ids);
 end
 
 perc_outside=[];
 perc_inside=[];
 
 for  k=1:length(data_matrix)
     
     column=data_matrix(:,k);
     column=column(~isnan(column));
     
     ix=find(column>limit_high);
     l=column(ix);
    ixx=find(column<limit_low);
     ll=column(ixx);
    outside=[ l;ll];

   p= (length(outside)/length(column))*100  ;   
    
   perc_outside=[perc_outside, p];
     
 end
 
 perc_outside=perc_outside(~isnan(perc_outside));
     
perc_outside_outside=[ perc_outside(2:5) perc_outside(16:end-1)]

perc_outside_inside=[perc_outside(6:15)]

perc_outside_pre=[perc_outside(1)]
perc_outside_post=[perc_outside(end)]

ranksum([perc_outside_pre perc_outside_post] , perc_outside_outside)
ranksum( perc_outside_inside, perc_outside_outside)
%% calculate Perc % outside for PBs inside and outside of Shuffled distributions

no_nan_pre_bird=pre_bird(~isnan(pre_bird));

perc_outside_fakes=[];

pb_category_outside=[3 4 5 6 17 18 20 23 28 33 38];
pb_category_inside=[7 8 9 10 11 12 13 14 15 16];

for q=1:1000
    
    perc_outside_fake=[];
 for  k=1:length(pb_category_outside)
     
     column=data_matrix(:,pb_category_outside(k));
     column=column(~isnan(column));
    
     fake_resp=[];
     
     for t=1:length(column)
         
     entry = no_nan_pre_bird(randi(length(no_nan_pre_bird)));
     
     fake_resp=[fake_resp entry];
     
     end
     
      ix=find(fake_resp>limit_high);
     l=fake_resp(ix);
    ixx=find(fake_resp<limit_low);
     ll=fake_resp(ixx);
    outside=[ l ll];

   p= (length(outside)/length(fake_resp))*100  ;   
    
   perc_outside_fake=[perc_outside_fake, p];
     
 end

     perc_outside_fakes=[perc_outside_fakes perc_outside_fake'];
    
end

fake_means=mean(perc_outside_fakes,2);

 signrank(mean(perc_outside_fakes,2), perc_outside_outside)
        
       
figure
plot([2*ones(1,length(perc_outside_outside))], (perc_outside_outside),'ro')
hold on
plot([2.5*ones(1,length(fake_means))], (fake_means),'bo')
hold on
plot([2*ones(1,length(fake_means));2.5* ones(1,length(fake_means))],[perc_outside_outside', fake_means]','-', 'Color', [0.5 0.5 0.5])
hold on
plot(0.25,perc_outside_pre , 'o')
hold on
plot(0.5,perc_outside_post , 'o')
hold on
plot(1,perc_outside_inside , 'ro')
xlim([0 3])
ylim([0 70])
% hold on
% plot([3*ones(1,length(fake_means))], (fake_means),'o')
% 

ylabel('% of whistle syllables outside central pitch range')
box off
xticks([0.25 0.5 1 1.5 2 2.5 ])
xticklabels({'Pre','Post', 'PBs inside','Shuffled PBs inside','PBs outside','Shuffled PBs outside'})
axis square
set(gca,'TickDir','out')

perc_inside_fakes=[];


for q=1:1000
    
    perc_inside_fake=[];
 for  k=1:length(pb_category_inside)
     
     column=data_matrix(:,pb_category_inside(k));
     column=column(~isnan(column));
    
     fake_resp=[];
     
     for t=1:length(column)
         
     entry = no_nan_pre_bird(randi(length(no_nan_pre_bird)));
     
     fake_resp=[fake_resp entry];
     
     end
     
      ix=find(fake_resp>limit_high);
     l=fake_resp(ix);
    ixx=find(fake_resp<limit_low);
     ll=fake_resp(ixx);
    outside=[ l ll];
    

   p= (length(outside)/length(fake_resp))*100  ;   
    
   perc_inside_fake=[perc_inside_fake, p];
     

 end

     perc_inside_fakes=[perc_inside_fakes perc_inside_fake'];
    
end

fake_means=mean(perc_inside_fakes,2);

 hold on
 plot([1.5*ones(1,length(fake_means))], (fake_means),'bo')
plot([1*ones(1,length(fake_means));1.5* ones(1,length(fake_means))],[perc_outside_inside', fake_means]','-', 'Color', [0.5 0.5 0.5])

 signrank(mean(perc_inside_fakes,2), perc_outside_inside)
 signrank(mean(perc_outside_fakes,2), perc_outside_outside)
 
 median(perc_outside_outside)
 
 ranksum(mean(perc_outside_fakes,2), mean(perc_inside_fakes,2))
 ranksum(perc_outside_inside, perc_outside_outside)
 ranksum(perc_outside_pre, mean(perc_inside_fakes,2))
 ranksum(perc_outside_pre, mean(perc_outside_fakes,2))
ranksum(perc_outside_post, mean(perc_inside_fakes,2))
 ranksum(perc_outside_post, mean(perc_outside_fakes,2))

 
 median(mean(perc_outside_fakes,2))
 
 
 %%   bin responses to each PB and calculate Perc %
 
 jj = [600:200:8600];
 
 PERC_pb_BIN=[];
 PERC_pre_BIN=[];
 range=100
 
 for t=7:16
     
 Y=data_matrix(:,(t));
 %Y = discretize(data_matrix(:,(t)),[900:200:9000]);
 
 Y=Y(~isnan(Y));
 
% figure
% histogram(Y,41)
 
 PERC_pb_BIN=[ PERC_pb_BIN, length(find(Y>=jj(t)-range & Y<=jj(t)+range))/length(Y)*100 ]

 W=data_matrix(:,1);
 W=W(~isnan(W));

PERC_pre_BIN=[PERC_pre_BIN, length( find(W>=jj(t)-range & W<=jj(t)+range))/length(W)*100 ]
 
 end
 
 signrank(PERC_pb_BIN, PERC_pre_BIN)

% PERC_pb_BIN - PERC_pre_BIN
 
 mean(PERC_pb_BIN)
 std(PERC_pb_BIN)
 median(PERC_pb_BIN)

 
% mean(PERC_pre_BIN)
% std(PERC_pre_BIN)
% median(PERC_pre_BIN)
oo=PERC_pb_BIN;
 o=PERC_pre_BIN;
 
figure(2222)
subplot(2,1,1)
 plot([2*ones(1,length(PERC_pb_BIN))], (PERC_pb_BIN),'o')
hold on
 plot([1*ones(1,length(PERC_pre_BIN))], (PERC_pre_BIN),'o')
plot([2*ones(1,length(PERC_pb_BIN));1* ones(1,length(PERC_pre_BIN))],[PERC_pb_BIN ;(PERC_pre_BIN)],'-', 'Color', [0.5 0.5 0.5])
xlim([0 3])
ylim([0 30])
axis square
box off
set(gca,'TickDir','out')
ylabel('% of accurate responses(Wsyll inside BINpb+-100Hz)')
box off
xticks([ 1 2  ])
xticklabels({'Pre','Responses'})


%%

 
 jj = [600:200:8600];
 
 PERC_pb_BIN=[];
 PERC_pre_BIN=[];
 range=100
%  7 8 9 10 11 12 13 14 15 16
 s=[3 4 5 6  17 18 20 23 28 33 38]
 
 for t=s
     
 Y=data_matrix(:,(t));
 %Y = discretize(data_matrix(:,(t)),[900:200:9000]);
 
 Y=Y(~isnan(Y));
 
% figure
% histogram(Y,41)
 
 PERC_pb_BIN=[ PERC_pb_BIN, length(find(Y>=jj(t)-range & Y<=jj(t)+range))/length(Y)*100 ]

 W=data_matrix(:,1);
 W=W(~isnan(W));

PERC_pre_BIN=[PERC_pre_BIN, length( find(W>=jj(t)-range & W<=jj(t)+range))/length(W)*100 ]
 
 end
 
 signrank(PERC_pb_BIN, PERC_pre_BIN)
 
  ranksum(PERC_pb_BIN, oo)


% PERC_pb_BIN - PERC_pre_BIN
 
% mean(PERC_pb_BIN)
% std(PERC_pb_BIN)
 median(PERC_pb_BIN)

 
% mean(PERC_pre_BIN)
% std(PERC_pre_BIN)
 median(PERC_pre_BIN)

 
 
figure(2222)
subplot(2,1,2)
 plot([2*ones(1,length(PERC_pb_BIN))], (PERC_pb_BIN),'o')
hold on
 plot([1*ones(1,length(PERC_pre_BIN))], (PERC_pre_BIN),'o')
plot([2*ones(1,length(PERC_pb_BIN));1* ones(1,length(PERC_pre_BIN))],[PERC_pb_BIN ;(PERC_pre_BIN)],'-', 'Color', [0.5 0.5 0.5])
xlim([0 3])
ylim([0 30])
axis square
box off
set(gca,'TickDir','out')
ylabel('% of accurate responses (Wsyll inside BINpb +- 100Hz)')
box off
xticks([ 1 2  ])
xticklabels({'Pre','Responses'})


%%
figure
 plot([1*ones(1,length(PERC_pb_BIN))], (PERC_pb_BIN),'o')
hold on
 plot([0.5*ones(1,length(PERC_pre_BIN))], (PERC_pre_BIN),'o')
plot([1*ones(1,length(PERC_pb_BIN));0.5* ones(1,length(PERC_pre_BIN))],[PERC_pb_BIN ;(PERC_pre_BIN)],'-', 'Color', [0.5 0.5 0.5])
 plot([1.5*ones(1,length(oo))], (oo),'o')
hold on
 plot([2*ones(1,length(o))], (o),'o')
plot([1.5*ones(1,length(o));2* ones(1,length(o))],[oo ;(o)],'-', 'Color', [0.5 0.5 0.5])

xlim([0 2.5])
ylim([0 35])
axis square
box off
set(gca,'TickDir','out')
ylabel('% of accurate responses (Wsyll inside BINpb +- 100Hz)')
xticks([ 0.5 1 1.5 2  ])
xticklabels({'Pre', 'Responses', 'Responses','Pre'})

%%
[p,h]=signrank([PERC_pb_BIN oo], [PERC_pre_BIN o])


std([PERC_pb_BIN oo])

median([PERC_pre_BIN o])

mean(PERC_pb_BIN)
mean(oo)

mean(PERC_pre_BIN)
mean(o)

%%
ff=[s 7:16];

figure
scatter([s 7:16], [PERC_pb_BIN oo])
hold on
scatter([s 7:16]+0.01, 1*[PERC_pre_BIN o])
xlim([1 41])
ylim([0 30])

yline(median([PERC_pb_BIN oo]), 'b--')
yline(median(1*[PERC_pre_BIN o]), 'r--')

ghg=[PERC_pb_BIN oo];
hgh=[PERC_pre_BIN o];

for k=1:length(ghg)
plot([ff(k)*ones(1,1);((ff(k)+0.01)* ones(1,1))],[ghg(k) ; hgh(k)],'-', 'Color', [0.5 0.5 0.5])

end

box off
set(gca,'TickDir','out')
ylabel('% of accurate responses (Wsyll inside BINpb +- 100Hz)')
xticks([sort([s 7:16] )] )
xticklabels({'1.0', '1.2', '1.4','1.6', '1.8', '2.0','2.2', '2.4', '2.6', '2.8','3.0',...
    '3.2', '3.4','3.6', '3.8','4.0', '4.4', '5.0','6.0', '7.0', '8.0'})

%%
q=[0 0.6 0.3 0 0.1 0.1 0.1 0 0.4 0.8 0.8 0.6 0.3 0.1 0 0.1 0.8 0.1 0.1 0 0]';
w=[4.9 0 4.1 3.8 2.7 3.3 2.4 1.9 1.2 1.4 2.4 1.6 2.2 3.2 2.3 1.5 1.9 2.2 2.7 1.7 1.5]';
e=[1.5 0.7 0 0.9 0.4 1.4 1.2 0.4 0.9 0.8 0.8 0.5 0.7 1.0 0.6 0.8 0.8 1.1 0.3 0.5 1.1]';
r=[7.6 7.4 10.0 0 9.6 7.9 5.2 6.2 5.9 2.8 10.2 4.7 10.1 6.0 5.9 4.0 6.4 1.9 3.7 11.5 2.1]';
t=[1.9 3.0 3.0 2.0 0 3.7 4.1 4.0 1.2 2.1 2.7 2.7 3.2 2.0 4.7 3.9 3.1 0.9 1.5 3.3 4.2]';
y=[7.8 10.1 8.7 11.8 10.9 0 7.7 9.3 4.3 6.1 7.1 8.2 6.5 6.6 11.8 10.4 7.4 6.2 7.2 12.2 10.4]';
u=[6.0 2.8 7.1 7.6 2.9 3.1 0 10.5 6.5 3.2 7.4 4.3 2.5 5.4 2.4 3.5 5.1 4.3 6.3 3.8 3.4]';
i=[1.8 2.0 2.4 2.2 1.9 1.5 3.7 0 3.7 2.2 2.2 0.9 0.8 2.5 2.3 2.1 3.1 1.8 1.3 1.4 2.1]';
o=[3.4 2.1 2.8 1.5 3.4 2.7 2.6 5.8 0 5.8 3.7 1.9 3.9 0.5 0.8 2.9 3.6 2.0 1.7 1.3 2.3]';
p=[23.0 16.2 12.5 10.1 13.4 10.7 6.8 16.5 28.5 0 15.8 16.3 17.4 12.0 3.8 8.0 12.1 15.5 4.2 10.0 11.5]';
a=[2.8 2.0 1.6 1.8 2.4 2.3 1.9 1.6 1.4 2.6 0 5.1 4.5 3.1 3.5 1.1 2.6 2.4 1.9 2.8 4.4]';
s=[9.4 9.1 6.5 8.0 11.4 7.1 7.3 7.2 6.2 9 15.2 0 13.1 11.7 12.7 8.3 8.6 8.7 6.8 11.4 13.8]';
d=[1.0 0.2 1.7 0.2 0.1 0.5 0 0.7 2.3 0.9 0.1 0.3 0 0.1 0 0 0.6 0.1 0.1 0.3 0.8]';
f=[4.1 4.2 4.3 3.8 4.6 5.1 3.9 2.4 2.0 3.4 3.6 4.2 3.5 0 9.7 5.2 5.4 4.2 4.0 5.0 6.0]';
g=[1.3 3.0 4.1 2.7 6.8 4.7 3.7 2.9 0.4 2.4 4.6 4.4 4.5 7.1 0 5.1 2.8 3.3 3.0 2.6 6.2]';
h=[6.0 5.3 7.3 4.3 4.9 6.6 9.6 8.6 6.9 5.1 4.1 6.6 3.2 6.1 8.5 0 6.9 4.3 4.5 6.7 6.3]';
j=[4.2 7.2 10.4 4.9 4.4 3.3 7.2 5.0 4.0 3.0 3.7 4.5 5.7 4.4 3.3 7.5 0 6.0 2.2 2.5 6.5]';
k=[6.5 8.0 6.3 4.8 4.3 5.0 7.2 5.9 7.6 3.6 2.0 2.5 6.7 4.2 5.9 6.4 8.1 0 9.8 3.9 7.0]';
l=[4.1 5.8 2.7 4.4 3.2 3.3 1.9 3.7 5.6 2.3 1.4 3.5 6.5 5.7 1.1 4.9 3.0 8.0 0 13.9 4.9]';
z=[1.8 1.2 1.6 3.8 2.5 3.2 1.9 2.1 0.9 1.2 2.2 3.0 2.0 1.6 3.3 1.6 1.4 1.1 1.7 0 1.9]';
x=[0.5 1.0 1.1 0.8 2.9 1.8 3.0 1.3 1.7 2.2 3.3 2.0 1.0 3.2 2.7 3.1 1.1 0.1 1.2 1.0 0]';

FNR_real=[q,w,e,r,t,y,u,i,o,p,a,s,d,f,g,h,j,k,l,z,x];

TPR_real=zeros(21)

tpr_real=[0.6 8.1 1.4 20.8 7.1 22.6 18.5 4.0 8.2 39.1 7.0 22.1 1.5 13.4 15.0 19.7 15.3 25.7 35.9 4.2 3.7];

for qq=1:21
TPR_real(qq,qq)=tpr_real(qq)
end

TPR_real(isnan(TPR_real))=0;

real_matrix=FNR_real+TPR_real;

q=    [2.6 2.6 3.0 3.2 2.9 3.1 3.1 2.4 1.3 1.7 1.3 2.8 3.2 2.4 1.8 3.2 2.9 2.9 2.2 2.2 3.4]';
   w= [8.4 9.8 9.2 10.4 10.8 11.1 10.2 8.0 12.0 9.9 10.8 9.5 10.9 8.1 9.2 9.1 9.6 9.0 10.0 9.2 9.4]';
 e =  [1.9 2.1 3.2 3.8 1.8 2.0 1.2 2.8 2.2 1.6 1.3 3.0 1.7 2.9 1.8 2.8 2.0 2.7 2.2 2.0 2.1]';
r  =  [10.9 9.2 9.5 7.9 8.3 8.3 8.4 10.2 8.4 9.4 9.9 7.8 7.6 9.8 8.0 8.4 9.8 9.5 7.6 9.3 8.1]';
  t=  [3.1 3.6 3.8 3.2 2.8 2.7 3.3 4.3 3.2 3.7 2.9 2.3 2.9 3.1 2.6 2.8 3.0 2.8 4.1 3.0 2.6]';
y  = [4.5 4.8 4.6 4.1 4.6 3.8 4.8 4.3 4.6 5.3 4.7 3.9 3.7 4.4 5.3 5.1 6.1 3.5 5.8 5.1 4.7]';
u  =  [2.3 2.1 2.4 2.7 2.5 3.7 1.6 2.2 3.2 2.3 1.4 2.4 2.3 3.2 3.0 2.7 2.5 2.3 2.1 3.3 3.1]';
 i =  [1.9 2.1 1.7 2.4 2.1 1.5 2.0 1.8 2.6 2.6 3.2 2.0 2.2 2.4 3.2 3.3 2.3 1.4 3.0 1.3 2.1]';
  o = [3.1 3.8 2.2 1.8 2.4 2.6 2.0 2.2 3.5 2.9 3.4 2.9 3.2 3.2 1.8 2.9 2.3 3.2 2.6 2.0 3.1]';
 p  = [15.9 15.0 19.3 18.4 16.2 19.3 18.1 20.2 14.3 20.3 17.0 17.9 18.0 17.2 18.1 15.7 15.3 18.1 14.4 17.0 18.8]';
  a = [3.1 3.0 3.5 3.7 3.7 2.7 3.4 4.0 4.3 4.0 4.6 3.2 4.4 4.5 3.2 3.1 3.4 3.3 3.7 3.8 4.2]';
s=    [13.8 14.9 12.2 11.6 13.7 12.9 12.5 11.2 13.7 12.7 13.5 14.6 13.9 15.4 13.6 14.4 15.6 14.2 13.9 14.3 11.8]';
d =   [1.0 1.5 1.3 1.1 0.6 1.0 1.4 0.7 1.3 0.8 1.4 1.3 1.3 0.9 1.4 1.3 1.0 1.4 1.0 1.2 1.1]';
 f =  [4.9 3.8 4.7 3.0 4.9 4.3 3.5 3.5 4.0 3.2 4.1 3.9 3.4 3.9 4.7 3.7 3.6 4.1 4.1 4.7 3.2]';
 g  = [1.5 1.9 0.9 1.8 0.9 1.5 1.9 1.5 1.6 1.4 2.2 1.7 1.2 2.0 2.4 1.9 1.6 2.0 1.5 2.5 2.8]';
 h  = [1.8 2.0 1.7 2.3 1.3 1.4 1.6 2.2 2.4 0.9 1.8 2.4 2.3 1.7 2.0 1.6 1.9 1.6 2.2 2.1 1.5]';
 j  = [5.7 2.5 4.0 4.2 4.0 2.7 4.1 3.8 3.6 2.8 4.4 3.7 3.7 3.6 3.2 3.6 2.6 2.8 4.2 3.4 1.8]';
 k  = [2.1 2.5 2.2 2.8 3.1 4.0 3.5 1.9 3.5 2.8 2.5 2.2 3.2 1.5 2.1 2.9 2.6 1.0 3.9 2.2 2.4]';
 l  = [6.3 4.9 4.0 5.4 5.9 5.5 6.4 5.6 5.0 4.9 4.2 4.9 4.2 4.5 4.5 6.0 6.6 6.8 4.7 4.8 6.0]';
z  =  [4.1 5.6 4.4 4.8 5.9 4.7 5.4 5.0 3.6 4.7 4.3 4.9 4.4 3.4 5.9 4.1 3.5 4.9 4.6 5.5 5.2]';
x  =  [1.3 2.1 2.2 1.7 1.8 1.3 1.5 2.2 1.9 2.1 1.1 2.6 2.5 2.0 2.3 1.3 1.9 2.4 2.3 1.0 2.6]';
    

shuffled_matrix=[q,w,e,r,t,y,u,i,o,p,a,s,d,f,g,h,j,k,l,z,x];

tpr_shuff=[ 2.6 9.8 3.2 7.9 2.8 3.8 1.6 1.8 3.5 20.3 4.6 14.6 1.3 3.9 2.4 1.6 2.6 1.0 4.7 5.5 2.6 ];



mean(tpr_real)
mean(tpr_shuff)
%%

figure
hold on
plot(real_matrix(11,:))


figure
imagesc([q,w,e,r,t,y,u,i,o,p,a,s,d,f,g,h,j,k,l,z,x])
cmap=flipud(copper(2048));
colormap(cmap);
cmap(1,:)=[1,1,1];
colormap(cmap)
colorbar
box off
axis square
set(gca,'TickDir','out')

%%
figure
imagesc(real_matrix)
cmap=(summer(2048));;
colormap(cmap);
cmap(1,:)=[1,1,1];
colormap(cmap)
colorbar
box off
box off
axis square
set(gca,'TickDir','out')
set(gca,'color','none')
xticks([1:21])
yticks([1:21])
xticklabels({'1.0', '1.2', '1.4','1.6', '1.8', '2.0','2.2', '2.4', '2.6', '2.8','3.0',...
    '3.2', '3.4','3.6', '3.8','4.0', '4.4', '5.0','6.0', '7.0', '8.0'})
yticklabels({'1.0', '1.2', '1.4','1.6', '1.8', '2.0','2.2', '2.4', '2.6', '2.8','3.0',...
    '3.2', '3.4','3.6', '3.8','4.0', '4.4', '5.0','6.0', '7.0', '8.0'})

caxis([0 60])

figure
imagesc(shuffled_matrix)
cmap=(summer);;
colormap(cmap);
%cmap(1,:)=[1,1,1];
colormap(cmap)
colorbar
box off
box off
axis square
set(gca,'TickDir','out')
set(gca,'color','none')
xticks([1:21])
yticks([1:21])
xticklabels({'1.0', '1.2', '1.4','1.6', '1.8', '2.0','2.2', '2.4', '2.6', '2.8','3.0',...
    '3.2', '3.4','3.6', '3.8','4.0', '4.4', '5.0','6.0', '7.0', '8.0'})
yticklabels({'1.0', '1.2', '1.4','1.6', '1.8', '2.0','2.2', '2.4', '2.6', '2.8','3.0',...
    '3.2', '3.4','3.6', '3.8','4.0', '4.4', '5.0','6.0', '7.0', '8.0'})
caxis([0 40])

%%
mean_diag=[];

ww=[-20:1:20];

for qq=1:41

x1 = diag(real_matrix,ww(qq));
mean_diag=[mean_diag, mean(x1)];

end


figure
subplot(1,2,1)
plot(ww,mean_diag)
box off
box off
axis square
set(gca,'TickDir','out')
ylim([0 16])
yline(100/21)

mean_diag=[];

ww=[-20:1:20];

for qq=1:41

x1 = diag(shuffled_matrix,ww(qq));
mean_diag=[mean_diag, mean(x1)];

end


subplot(1,2,2)
plot(ww,mean_diag)
box off
box off
axis square
set(gca,'TickDir','out')
ylim([0 16])
yline(100/21)


%%

bird_stim_id=[pb_bird, pb_stim];

[q,qq]=find(pb_stim==2400 |  pb_stim==2800 |  pb_stim==3200 |  pb_stim==5000 |  pb_stim==6000)

limited_bird_pb=[pb_bird(q), pb_stim(q)]

%%

q=[24.2 26.0 31.4 11.1 7.4];
w=[11.6 50.0 26.4 6.1 6.0];
e=[9.5 22.2 53.9 6.7 7.7];
r=[9.3 22.7 27.2 32.2 8.6];
t=[14.0 10.0 22.2 14.1 39.7];

mm=[q;w;e;r;t]

figure
imagesc(mm)
cmap=(summer);;
colormap(cmap);
%cmap(1,:)=[1,1,1];
colormap(cmap)
colorbar
box off
box off
axis square
set(gca,'TickDir','out')
set(gca,'color','none')
caxis([0 60])

yticks([1:5])
xticklabels({'2400',  '2800',...
    '3200', '5000','6000'})
yticklabels({'2400',  '2800',...
    '3200', '5000','6000'})

mean_diag=[];

ww=[-20:1:20];

for qq=1:41

x1 = diag(mm,ww(qq));
mean_diag=[mean_diag, mean(x1)];

end


figure
plot(ww,mean_diag)
box off
box off
axis square
set(gca,'TickDir','out')


%%
x1 =[24.2000 50.0000 53.9000 32.2000 39.7000];
x2 =[80.6000 45.5000 77.8000 22.2000 92.3000];

median(x2)
median(x1)
[p,h]=signrank(x1,x2)