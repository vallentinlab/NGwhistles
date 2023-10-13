%% Whistle analysis 2.0 for PAIRS %%
clear all
close all
clc

% Read all bird's data...
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*ed.csv');
whistle_table = [];

%%
low_limit = 1.0e+03 *[1.1852,2.5794,6.1129,3.9992,1.7282,2.2298,3.2687,1.5235,5.2733,3.5808,7.4114,5.6278,1.9587,6.7073,4.4758,2.9400,8.1716,4.8752]

high_limit =1.0e+03 *[1.5231,2.9382,6.6943,4.4624,1.9568,2.5741,3.5783,1.7256,5.6194,3.9623,8.1415,6.0983,2.2199,7.3649,4.8667,3.2649,9.1690,5.2637]


centroids =1.0e+03 *[ 1.4339,1.6141,1.8399,2.0739,2.3733,2.7857,3.0956,3.4329,3.7155,4.2739,4.6662,5.0836,5.4598,5.7941,6.4262,6.9774,7.7834,8.5399]

low_limit= 1000*[5.0780  2.9580 1.1852  2.2022 6.1486 1.7882 4.0919 7.4114     2.5834 3.3563 ]  ;           
                  
high_limit = 1000*[ 6.1365  3.3522  1.7858 2.5794 7.3649 2.1965 5.0671 9.1690  2.9538  4.0763    ];    
    

centroids =1.0e+03 *[1.5515 5.5585 3.1223 8.0409 2.3710  2.7893  2.0246 6.7185  4.5893    3.5879];
    
    
    
    
   
   
    
   
 

all_bird_id=[];
bird_pitch=[];
matrix_freq = nan(600);
matrix_id=nan(600);
%% Load data/ adress folder %%
for bird=1:20
    
if bird == 1
    fprintf("Count with me! :)\n");
end
bird

filename= char(strcat(pathname,'\',filelist(bird,1).name));
[data,text] = xlsread(filename);

%% Read Excel %%

bird_freq = data(:,5);

bird_pitch=[bird_pitch;bird_freq];

all_bird_id=[all_bird_id;bird*ones(length(bird_freq),1)];


 matrix_freq(1:length(bird_freq),bird)=bird_freq;
 matrix_id(1:length(bird_freq),bird)=bird*ones(length(bird_freq),1);
 
figure(1)

  oo=[rand rand rand];
  scatter(bird_freq,bird*ones(length(bird_freq),1),...
        'MarkerEdgeColor', 'none','MarkerFaceColor',oo);
    for t=1:length(low_limit)
        xline(low_limit(t), '--')
    hold on
    end
     hold on
     
  

 
end % for loop, all birds


%%
 
 J=customcolormap([0 0.5 1], [0 0.2 0; 0 0.5 0;1 1 1 ]);

figure
hold on
 [values, centers] = hist3([ bird_pitch, all_bird_id],[250 40]);
imagesc(centers{:}, values.');
colorbar
axis xy 
xlim([1000 10000]);
    ylim([0 22]);
         for t=1:length(low_limit)
        xline(low_limit(t), '--')
    hold on
    end
 colormap(J)

box off
set(gca,'TickDir','out')

hold on
for ii=1:length(centroids)
    plot(centroids(ii),21, '*', 'MarkerSize', 10)
hold on
end


%%


numb=[10 9 10 10 9 10 9 10 10 10 10 10 9 10 10 10 10 10 10 10];



length(numb)

mean(numb)
std(numb)
median(numb)