clear all; close all; clc

[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);
filelist = dir('*.xlsx');

PRE_pitch=[]; PRE_set=[]; PRE_birdID=[]; PRE_snippet=[]; PRE_whistle_syll_ID =[];
POST_pitch=[]; POST_set=[]; POST_birdID=[]; POST_snippet=[]; POST_whistle_syll_ID=[];
PB_pitch_bird = []; PB_pitch_stim = []; PB_birdID=[]; PB_snippet=[]; PB_whistle_syll_ID =[];

for bird=1:11

bird

filename= char(strcat(pathname,'\',filelist(bird,1).name));
[~,sheet_name]=xlsfinfo(filename);
data=[];
for k=1:numel(sheet_name)
  [num,txt,raw]=xlsread(filename,sheet_name{k});
data{k}=raw;
dataa{k}=num;
raw=[];
num=[];
end

PRE_pitch=[PRE_pitch; dataa{1,1}(:,3)];
PRE_set=[ PRE_set; dataa{1,1}(:,10) ];
PRE_birdID=[PRE_birdID zeros*(1:length(dataa{1,1}(:,10)))+bird];
PRE_snippet=[PRE_snippet;dataa{1,1}(:,8) ];

B=convertStringsToChars(data{1,1}(:,4));
  newStr = erase(B(1:end,1),"w"); whistle_id=[];
  for w=1:length(newStr)
  C=str2num(newStr{w});
   whistle_id=[whistle_id; C];
  end
  
  PRE_whistle_syll_ID=[PRE_whistle_syll_ID; whistle_id];
  
  
PB_pitch_bird = [PB_pitch_bird; dataa{1,2}(:,3)     ];
PB_pitch_stim = [PB_pitch_stim; dataa{1,2}(:,14)      ];
PB_birdID=[PB_birdID zeros*(1:length(dataa{1,2}(:,10)))+bird     ];
PB_snippet=[PB_snippet; dataa{1,2}(:,8)      ];

B=convertStringsToChars(data{1,2}(:,4));
  newStr = erase(B(1:end,1),"w"); whistle_id=[];
  for w=1:length(newStr)
  C=str2num(newStr{w});
   whistle_id=[whistle_id; C];
  end
  PB_whistle_syll_ID=[PB_whistle_syll_ID; whistle_id];


POST_pitch=[POST_pitch; dataa{1,3}(:,3)];
POST_set=[ POST_set;dataa{1,3}(:,10) ];
POST_birdID=[POST_birdID zeros*(1:length(dataa{1,3}(:,10)))+bird];
POST_snippet=[POST_snippet;dataa{1,3}(:,8) ];

B=convertStringsToChars(data{1,3}(:,4));
  newStr = erase(B(1:end,1),"w"); whistle_id=[];
  for w=1:length(newStr)
  C=str2num(newStr{w});
   whistle_id=[whistle_id; C];
  end
  POST_whistle_syll_ID=[POST_whistle_syll_ID; whistle_id];
  
end

PRE_birdID=PRE_birdID' ; PB_birdID=PB_birdID'; POST_birdID=POST_birdID';
%%



clearvars -except PRE_pitch PRE_set PRE_birdID PRE_snippet PRE_whistle_syll_ID ...
POST_pitch POST_set POST_birdID POST_snippet POST_whistle_syll_ID ...
PB_pitch_bird  PB_pitch_stim  PB_birdID PB_snippet PB_whistle_syll_ID


%% autocorr whistles in the PRE - looking for history of whistle songs

% for q=1:11
%     s=[]; ss=[]; w=[]; ww=[]; id_pre=[]; idd_pre=[];
%     [s,ss]=find(PRE_birdID==q);
%     
%     B=(PRE_set(s));
%     
%     [w,ww]=find(B==2);
% 
%     A=(PRE_pitch(s));
%     
%     C=A(w);
%     
%     D=PRE_whistle_syll_ID(s);
%     
% [id_pre, idd_pre]=find(D(w)==1);
% 
% 
% n=length(C(id_pre));
% 
% % figure
% % C(id_pre)
% % autocorr(C(id_pre), n-1)
% 
%  [acf,lags] = autocorr(C(id_pre), n-1);
% figure(q)
% plot(acf,'r', 'LineWidth', 5)
% hold on
% YourVector=C(id_pre);
% for o=1:1000
% YourVector=YourVector(randperm(length(YourVector)));
% [acf,lags] = autocorr(YourVector, n-1);
% 
% plot(acf,'k')
% end
% 
% [acf,lags] = autocorr(C(id_pre), n-1);
% figure(q)
% plot(acf,'r', 'LineWidth', 5)
% 
% figure
% 
% hold on
% YourVector=C(id_pre);
% for o=1:1000
% YourVector=YourVector(randperm(length(YourVector)));
% [z,zz] = xcorr(YourVector, 'normalized');
% 
% plot(zz,z, 'k')
% 
% end
% [z,zz]=xcorr(C(id_pre),'normalized');
% plot(zz,z,'r', 'LineWidth', 5)
% end
% 
% figure
% plot(C(id_pre))
%%

%% History of Delta of whistles

original_p=[];
diff_p=[];

L=[];
for q=3%1:11
    s=[]; ss=[]; w=[]; ww=[]; id_pre=[]; idd_pre=[];
    
    [s,ss]=find(PRE_birdID==q);
    B=(PRE_set(s));
    [w,ww]=find(B==2);
    A=(PRE_pitch(s));
    C=A(w);
    D=PRE_whistle_syll_ID(s);
    
[id_pre, idd_pre]=find(D(w)==1);

Q=C(id_pre);

figure(1000)
% subplot(3,1,1)
% scatter((Q(1:end-1)),(Q(2:end)))
% [rh,pval_or] = corr((Q(1:end-1)),(Q(2:end)), 'Type','Pearson','Rows','complete')
% axis square
% xlabel('Pitch (i)')
% ylabel('Pitch (i+1)')
% original_p=[original_p, pval_or];

subplot(11,2,q+q-1)
scatter(diff(Q(1:end-1)),diff(Q(2:end)))
[rhoo,pval_diff] = corr(diff(Q(1:end-1)),diff(Q(2:end)), 'Type','Pearson','Rows','complete')
axis square
box off
title(strcat('bird #', num2str(q)))
 xlim([-6000 6000])
 ylim([-6000 6000])
xlabel('wp(n)-wp(n+1) (Hz)')
ylabel('wp(n+1)-wp(n+2) (Hz)')


diff_p=[diff_p, pval_diff];

F=[];
YourVector=((C(id_pre)));
for o=1:1000%:1000
YourVector=YourVector(randperm(length(YourVector)));

YY=diff(YourVector);

[rho,pval] = corr(YY(1:end-1),YY(2:end), 'Type','Pearson','Rows','complete');

F=[F, rho];
end


subplot(11,2,q+q)
histogram(F,1000,'Normalization','probability')
xline(rhoo, 'r', 'LineWidth', 3)

box off
title(strcat('bird #', num2str(q)))
 xlim([-1 1])
 ylim([0 0.02])
ylabel('Probability')
xlabel('Correlation value')


L=[ L (length(find(F<rhoo)))/1000];

end

L
%%

for q=1:11
    s=[]; ss=[]; w=[]; ww=[]; id_pre=[]; idd_pre=[];
    
    [s,ss]=find(PRE_birdID==q);
    B=(PRE_set(s));
    [w,ww]=find(B==2);
    A=(PRE_pitch(s));
    C=A(w);
    D=PRE_whistle_syll_ID(s);
    
[id_pre, idd_pre]=find(D(w)==1);

Q=C(id_pre);
figure (1000)
subplot(2,6,q)
autocorr((Q), 11);
axis square
% subplot(1,2,1)
% plot(Q)
axis square
box off
title(strcat('bird #', num2str(q)))
set(gca,'TickDir','out')
grid off

xlim([-1 12])
ylim([-1.5 1.5])
ylabel('Correlation')
xlabel('Whistle position prior to focal whistle')

end

%%

%% sollwert - first to last syllable of PREs

PRE=[PRE_pitch PRE_whistle_syll_ID PRE_set PRE_birdID];

[e,ee]=find(PRE(2:end,2)==1);

FIRSTs= PRE_pitch(e);
LASTs=PRE_pitch(e-1);

F=FIRSTs-FIRSTs;
L=LASTs-FIRSTs;

figure
subplot(1,3,1)
 plot([1*ones(1,length(F))], (F),'o')
hold on
 plot([2*ones(1,length(L))], (L),'o')
plot([1*ones(1,length(F));2* ones(1,length(L))],[F, (L)]','-', 'Color', [0.5 0.5 0.5])
xlim([0 3])
ylim([-500 500])
yline(0, 'k--')
axis square
box off
set(gca,'TickDir','out')
ylabel('Delta Pitch from first syllable (Hz)')
box off
xticks([ 1 2  ])
xticklabels({'First syll','Last syll'})

figure
scatterplot(L)


% figure
% scatter(PRE_pitch(e),(L) )
% ylim([-1000 1000])
% yline(0, 'k--')
% axis square
% box off
% set(gca,'TickDir','out')
% ylabel('Delta Pitch from first syllable (Hz)')
% xlabel('Pitch of first syll')
% [rho,pval] = corr(PRE_pitch(e),(L), 'Type','Pearson','Rows','complete')
% 


[r,rr]=find(L>0);
subplot(1,3,2)

scatter(FIRSTs(r),L(r) )
[rho,pval] = corr(FIRSTs(r),L(r), 'Type','Pearson','Rows','complete')
[r,rr]=find(L<0);
hold on
scatter(FIRSTs(r),L(r) , 'r')
[rho,pval] = corr(FIRSTs(r),L(r), 'Type','Pearson','Rows','complete')
yline(0, 'k--')
axis square
box off
set(gca,'TickDir','out')
ylabel('Delta Pitch from first syllable (Hz)')
xlabel('Pitch of first syll')
ylim([-500 500])

subplot(1,3,3)

scatter(PRE_pitch(e),abs(L) )
ylim([0 500])
yline(0, 'k--')
axis square
box off
set(gca,'TickDir','out')
ylabel('Delta Pitch from first syllable (Hz)')
xlabel('Pitch of first syll')
[rho,pval] = corr(PRE_pitch(e),abs(L), 'Type','Pearson','Rows','complete')
length(L)

%% real distirbution of responses to calculate deltas from PBs

jj = [600:200:8600];

data_matrix = nan(10000);
data_matrix(1:length(PRE_pitch),1) = PRE_pitch;
data_matrix(1:length(POST_pitch),41) = POST_pitch;
 for ww=2:length (jj)
  ids=find(PB_pitch_stim==jj(ww));
  data_matrix(1:length(ids),ww)=PB_pitch_bird(ids);
 end
 
 STDs=[];
 
 figure 
 second_peak=[];
peak_all_bird=[];
peak_all_pb=[];
% Playback stimuli
for qq=1:length(jj);
    subplot(1,length(jj),qq);
    ids=find(PB_pitch_stim==jj(qq));
    PB_pitch_bird(ids);
    if (sum(PB_pitch_stim(ids))>0)==1;
        s2 = scatter_kde(PB_pitch_stim(ids)+(10+(rand(length(ids),1))*30),PB_pitch_bird(ids));
        s2.SizeData = 15;
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(jj(qq)));
        yline(mean(PB_pitch_stim(ids)),'-','LineWidth',7); % pb frequency
      
        [gg,bin]=hist(PB_pitch_bird(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
        hold on
        [ff, Xii, u] = ksdensity(data_matrix(:,qq),'Support','positive');
        [a,b]=findpeaks(ff);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xii(cc(end,2));
        yline(freq1,'-','LineWidth',7,'Color', 'r' ) % first peak
        
        if (length(a)>1)==1;
          freq2=Xii(cc(end-1,2));
              
          %yline(freq2,'-','LineWidth',8,'Color', 'm' );
          
        else
            freq2=nan;
        end
peak_all_pb=[peak_all_pb;jj(qq)];
        peak_all_bird=[peak_all_bird;freq1];
        
        second_peak=[second_peak; freq2];
    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(jj(qq)));
        hold on
    end
hold on


STDs=[STDs; nanstd(PB_pitch_bird(ids))];

end
STDs

%% compare STDs Gambia Germany

STDs=STDs(~isnan(STDs));

STDs([8,10,12,18,19])
gambia=[643.312358316530;419.422111337561;170.172744900527;1050.64364267983;540.529893187615];


median(STDs)
median(gambia)

[p,h]=ranksum(STDs, gambia)

%% octave idea


% octt=[peak_all_bird./peak_all_pb peak_all_pb./peak_all_bird]
% oct=[];
% 
% for w=1:length(octt)
%     
% oct=[oct max(octt(w,:))]
% 
% end
% 
% octttt=[peak_all_bird_shuf./peak_all_pb peak_all_pb./peak_all_bird_shuf]
% octtt=[];
% 
% for w=1:length(octttt)
%     
% octtt=[octtt max(octttt(w,:))]
% 
% end
% 
% figure;
% plot(peak_all_pb,oct,'o')
% 
% hold on
% 
% plot(peak_all_pb,octtt, 'or')
% 
% yline(1, '--k')
% yline(2, '--k')
% yline(3,'--k')
% ylim([0 4])


%% now shuffle!

data_matrix_shuff = nan(10000);
data_matrix_shuff(1:length(PRE_pitch),1) = PRE_pitch;
data_matrix(1:length(POST_pitch),41) = POST_pitch;

PB_pitch_bird_shuff=PB_pitch_bird(randperm(length(PB_pitch_bird)));
% 

[p,h]=ranksum(PB_pitch_bird_shuff, PRE_pitch)

%  for ww=2:length (jj)
%      
%   ids=find(PB_pitch_stim==jj(ww));
%  data_matrix_shuff(1:length(ids),ww)=PB_pitch_bird(randperm(length(PB_pitch_bird))); %%random sampling from the shuffled responses
%   
%   
%     % data_matrix_shuff(1:length(ids),ww)=PRE_pitch(randperm(length(ids))); %% random sampling from the PRE
% 
%  end

figure
peak_all_bird_shuf=[];
peak_all_pb_shuf=[];
% Playback stimuli
subplot (1,length(jj),1);
ylim([0 9000]);

for qq=2:length(jj);
    subplot(1,length(jj),qq);
    ids=find(PB_pitch_stim==jj(qq));
    PB_pitch_bird_shuff(ids);
    if (sum(PB_pitch_stim(ids))>0)==1;
        s2 = scatter_kde(PB_pitch_stim(ids)+(10+(rand(length(ids),1))*30),PB_pitch_bird_shuff(ids));
        s2.SizeData = 15;
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(jj(qq)));
        yline(mean(PB_pitch_stim(ids)),'-','LineWidth',7); % pb frequency
      
        [gg,bin]=hist(PB_pitch_bird_shuff(ids),25);
        [r_max,b]=max(gg);
        pb_max=bin(b);
        hold on
        [f, Xi, u] = ksdensity(PB_pitch_bird_shuff(ids),'Support','positive');
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        yline(freq1,'-','LineWidth',7,'Color', 'r' ) % first peak
peak_all_pb_shuf=[peak_all_pb_shuf;jj(qq)];
        peak_all_bird_shuf=[peak_all_bird_shuf;freq1];
    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(jj(qq)));
        hold on
    end
hold on
end

%%



for ww=2:length (jj)
     
  ids=find(PB_pitch_stim==jj(ww));
     data_matrix_shuff(1:length(ids),ww)=PRE_pitch(randperm(length(ids))); %% random sampling from the PRE

end



figure
 peak_all_bird_shuf=[];
peak_all_pb_shuf=[];
for qq=1:length(jj);
    subplot(1,length(jj),qq);
    ids=find(PB_pitch_stim==jj(qq));
    PB_pitch_bird_shuff=PRE_pitch(randperm(length(ids)));
    
hold on
    if (sum(PB_pitch_stim(ids))>0)==1;
        s2 = scatter_kde(PB_pitch_stim(ids)+(10+(rand(length(ids),1))*30),PB_pitch_bird_shuff);
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(jj(qq)));
        yline(mean(PB_pitch_stim(ids)),'-','LineWidth',7); % pb frequency
      
%         [gg,bin]=hist(PB_pitch_bird_shuff(ids),25);
%         [r_max,b]=max(gg);
%         pb_max=bin(b);
        hold on
        [f, Xi, u] = ksdensity(PB_pitch_bird_shuff,'Support','positive');
        [a,b]=findpeaks(f);
        c=[a;b];
        cc=sortrows(c');
        freq1=Xi(cc(end,2));
        yline(freq1,'-','LineWidth',7,'Color', 'r' ) % first peak
peak_all_pb_shuf=[peak_all_pb_shuf;jj(qq)];
        peak_all_bird_shuf=[peak_all_bird_shuf;freq1];
    else 
        ylim([0 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(jj(qq)));
        hold on
    end
hold on
end
 



%%



figure
for qq=1:5%length(jj);
     if (nansum(data_matrix(:,qq))>0)==1;
[ff, Xii, u] = ksdensity(data_matrix(:,qq),'Support','positive');

hold on
plot(Xii,ff)
     end
end
 
plot(Xii,(ff-f))
yline(0)


1+nan



%% plot and test!

delta=peak_all_bird-peak_all_pb;
delta_shuff=peak_all_bird_shuf-peak_all_pb_shuf;

figure
subplot(2,2,1)
stem(peak_all_pb,delta)
xlim([0 9000])
ylim([-6000 3000])
xline(1700,'--')
xline(3700,'--')
box off
axis square
set(gca,'TickDir','out')
% yticks([-1000 -500 0 500 1000])
ylim([-6500 3000])
ylabel('Delta pitch (Max-Pb) (Hz)')
xlabel('Pitch of whistle playback (Hz)')

figure

subplot(1,3,1)
stem(peak_all_pb_shuf,delta_shuff, 'r')
xlim([0 9000])
ylim([-6000 3000])
xline(1700,'--')
xline(3700,'--')
box off
axis square
set(gca,'TickDir','out')
% yticks([-1000 -500 0 500 1000])
ylim([-6500 3000])
ylabel('Delta pitch (Max-Pb) (Hz)')
xlabel('Pitch of whistle playback (Hz)')

median(abs(delta_shuff(5:14)))

% signrank(abs(delta(5:14)), abs(delta_shuff(5:14)))
% signrank((abs([delta(1:4) ;delta(15:end)])), (abs([delta_shuff(1:4) ;delta_shuff(15:end)]))

subplot(1,3,2)
 plot([1*ones(1,length(delta(5:14)))], abs(delta(5:14)),'o')
hold on
 plot([2*ones(1,length(delta_shuff(5:14)))], abs(delta_shuff(5:14)),'ro')
plot([1*ones(1,length(delta_shuff(5:14)));2* ones(1,length(delta_shuff(5:14)))],...
    [abs(delta(5:14))' ;abs(delta_shuff(5:14))'],'-', 'Color', [0.5 0.5 0.5])
xlim([0 3])
ylim([0 3000])
axis square
box off
yline(0, 'k--')
set(gca,'TickDir','out')
ylabel('Delta pitch (Max-Pb) for PBs inside central range (Hz)')
box off
xticks([ 1 2  ])
xticklabels({'Responses','Shuffled'})

subplot(1,3,3)
 plot([1*ones(1,length([delta(1:4) ;delta(15:end)]))], ([abs(delta(1:4)) ;abs(delta(15:end))]),'o')
hold on
 plot([2*ones(1,length([delta_shuff(1:4) ;delta_shuff(15:end)]))], ([abs(delta_shuff(1:4)) ;abs(delta_shuff(15:end))]),'ro')
plot([1*ones(1,length([delta_shuff(1:4) ;delta_shuff(15:end)]));2* ones(1,length([delta_shuff(1:4) ;delta_shuff(15:end)]))],...
    [[abs(delta(1:4)) ;abs(delta(15:end))]' ;([abs(delta_shuff(1:4)) ;abs(delta_shuff(15:end))])'],'-', 'Color', [0.5 0.5 0.5])
xlim([0 3])
ylim([0 3000])
yline(0, 'k--')
axis square
box off
set(gca,'TickDir','out')
ylabel('Delta pitch (Max-Pb) for PBs outside central range (Hz)')
box off
xticks([ 1 2  ])
xticklabels({'Responses','Shuffled'})


%%

figure
scatter(peak_all_pb,delta)
hold on
scatter(peak_all_pb_shuf,delta_shuff, 'r')

%%


%%

wwx=[3600:200:9000];

data_matrix=nan(10000);


for w=2:length (wwx)
    ids=find(PB_pitch_stim==wwx(w));
    data_matrix(1:length(ids),w)=PB_pitch_bird(ids);
 
end

whistle_over_pref=[];
stim_over_pref=[];
peak_high_only=[];
figure
for qq=1
    subplot(1,length(wwx)+1,qq)
   
    xlim([-10 50]);
        ylim([3652 9000]);
   
end
for qq=2:length(wwx);
    subplot(1,length(wwx)+1,qq);
    ids=find(PB_pitch_stim==wwx(qq));
    ddd=PB_pitch_bird(ids);
    ix=find(ddd>3652);
     l=ddd(ix);
    lll=[l];
    if (sum(PB_pitch_stim(ids))>0)==1;
        scatter_kde(PB_pitch_stim(ids(1:length(lll)))+(10+(rand(length(ids(1:length(lll))),1))*30),lll);
        ylim([3652 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwx(qq)));
        yline(mean(PB_pitch_stim(ids)),'-','LineWidth',8);
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
               peak_high_only=[peak_high_only;freq1];

    else 
        ylim([3652 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwx(qq)));
        hold on
    end
    hold on
    stim_over_pref=[stim_over_pref; PB_pitch_stim(ids(1:length(lll)))];
    whistle_over_pref=[whistle_over_pref; lll];
end


%%

  
    %%
    

    
   stim_over_pref_shuff=stim_over_pref(randperm(length(stim_over_pref)));
    whistle_over_pref_shuff= whistle_over_pref(randperm(length( whistle_over_pref)));

    figure
    scatter(stim_over_pref_shuff,whistle_over_pref_shuff)
    
    figure
for qq=1
    subplot(1,length(wwx)+1,qq)
   
    xlim([-10 50]);
        ylim([3652 9000]);
   
end

peak_high_only_shuff=[];
for qq=2:length(wwx);
    subplot(1,length(wwx)+1,qq);
    ids=find(stim_over_pref_shuff==wwx(qq));
    ddd=whistle_over_pref_shuff(ids);
    ix=find(ddd>3652);
     l=ddd(ix);
    lll=[l];
    if (sum(stim_over_pref_shuff(ids))>0)==1;
        scatter_kde(stim_over_pref_shuff(ids(1:length(lll)))+(10+(rand(length(ids(1:length(lll))),1))*30),lll);
        ylim([3652 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwx(qq)));
        yline(mean(stim_over_pref_shuff(ids)),'-','LineWidth',8);
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
               peak_high_only_shuff=[peak_high_only_shuff;freq1];

    else 
        ylim([3652 9000]);
        set(gca,'xtick',[],'ytick',[]);
        xlabel(convertCharsToStrings(wwx(qq)));
        hold on
    end
    hold on
   
end


%%
pbck=[38 40 44 50 60 70 80]*100;

delta_shuff=pbck-peak_high_only_shuff'

delta_=pbck-peak_high_only'

figure(22222)
subplot(1,2,1)
scatter(pbck,delta_shuff,'b')
 hold on
yline(0)
scatter(pbck,delta_,'r')
axis square
set(gca,'TickDir','out')
ylabel('Delta pitch (Max-Pb) for PBs outside central range (Hz)')
xlabel('Pitch of whistle playback (Hz)')


[h,p]=signrank(abs(delta_shuff(1:5)), abs(delta_(1:5)))
[h,p]=signrank(abs(delta_shuff), abs(delta_))

[h,p]=ttest((delta_shuff(1:5)) )


figure(22222)
subplot(1,2,2)

 plot([1*ones(1,length([delta_]))], ([abs(delta_)]),'ro')
hold on
 plot([2*ones(1,length([delta_shuff]))], ([abs(delta_shuff)]),'bo')
plot([1*ones(1,length([delta_shuff]));2* ones(1,length([delta_shuff]))],...
    [abs(delta_); abs(delta_shuff)],'-', 'Color', [0.5 0.5 0.5])
xlim([0 3])
ylim([0 5000])

yline(0, 'k--')
axis square
box off
set(gca,'TickDir','out')
ylabel('Delta pitch (Max-Pb) for PBs outside central range (Hz)')
box off
xticks([ 1 2  ])
xticklabels({'Observed','Shuffled'})
set(gca,'TickDir','out')
