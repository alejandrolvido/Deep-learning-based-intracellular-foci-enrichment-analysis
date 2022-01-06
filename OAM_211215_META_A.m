



path_h='C:\Users\aleja\OneDrive\Documents\MATLAB\MATLAB_doc\Scores_Home\';
folr='Vac3_Bck\';
exp_name='211124_cdc14_atg13_foci';
METAD=cell(9,1);

for pos=1:405

load([path_h folr exp_name '_pos' num2str(pos) '_extr_gradient'])

Gnumber=9;

for itl=1:length(cell_Total_mNG)

   if    ~isempty((cell_suma_mNG{itl}))
    
   gfocci=find(~isnan(cell_suma_mNG{itl}(:,1)));
     metadat=cell(9,1);
    for iu=1:length(gfocci)
    metadat{Gnumber}(iu,1)=pos;
    metadat{Gnumber}(iu,2)=itl;
    metadat{Gnumber}(iu,3)=gfocci(iu);
    metadat{Gnumber}(iu,4)= cell_Total_mNG{itl};
    metadat{Gnumber}(iu,5)=cell_Total_Ru{itl};
    
    metadat{1}(iu,:)=cell_suma_mNG{itl}(iu,:);
    metadat{2}(iu,:)=cell_averag_mNG{itl}(iu,:);
    metadat{3}(iu,:)=cell_maxima_mNG{itl}(iu,:);
    metadat{4}(iu,:)=cell_variation_mNG{itl}(iu,:);
    metadat{5}(iu,:)=cell_suma_Ru{itl}(iu,:);
    metadat{6}(iu,:)=cell_averag_Ru{itl}(iu,:);
    metadat{7}(iu,:)=cell_maxima_Ru{itl}(iu,:);
    metadat{8}(iu,:)=cell_variation_Ru{itl}(iu,:);
  
    end
    
    
%      for iyw=1:length(metadat)
      
    for ixr=1:Gnumber  
        if pos==1 && itl==1
        METAD{ixr}=metadat{ixr};
        else
        METAD{ixr}=[METAD{ixr}; metadat{ixr}];
        end
    end

%      end
   end
                      
                      
end

end

save('METAD','METAD')

for iuy=[2 6]%1:Gnumber-1
figure 
imagesc(METAD{iuy});colorbar%;caxis([-5000 5000])
% figure
% plot(nanmean(METAD{iuy}))
% title(num2str(iuy))
end

ecells=find(METAD{9}(:,1)<=79);
fcells=find(METAD{9}(:,1)>79);

METAD{9}(ecells,6)=1;
METAD{9}(fcells,6)=2;

for iy=1:8
    figure
for ih=1:2
  bcells=find(METAD{9}(:,6)==ih);
plot(nanmean(METAD{iy}(bcells,:)));
hold on
end
hold off
end


kgroup=4;
pooledMat_noNaN =METAD{2}(:,1:2);
[grId,C,sumd] = kmeans(pooledMat_noNaN, kgroup, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance','cityblock');
tabulate(grId)
METAD{9}(:,7)=grId;

for iw=1:2
METAD2=cell(size(METAD));
negpos=iw;
for iig=1:9
    
METAD2{iig}=METAD{iig}((METAD{9}(:,6)==negpos),:);

end

% define in previous step
SPL1=METAD2;
colrs=['r'; 'b'; 'k'; 'g'; 'c';'m';'y';'b'; 'r'; 'k'; 'g'; 'c';'m';'y']; % for plotting
limX=([1 20]); %endti
limY=57;
group=(1:kgroup);% specify number of clusters to plot 
Kgroup=kgroup;
axis=[0 200]; % adjust according to marker intensity 
AAA=0;
BBB=20-AAA;%endti;%25;
sr=1; 
sr1=1;
Snumber=0;
Snumber2=0;
nn2=1;
CL_ID=7;
starti=0; % first time point80
endti=20; % last time point
Gnumber=9;
sae=0;

for ii=[6]
      % ii=2
%%% asign marker to analyse as Cdc10-arranged matrix
    timeseriesM=SPL1{ii}; 
    titel1='Intensity ';char(SPL1{ii});  
%%            
% calculate mean curves using each lane of the microfluidic device as
% biological replicate
numCells = 1:length(SPL1{ii}(:,nn2));
chamId = SPL1{Gnumber}(:,nn2); 
CH_no=max(SPL1{Gnumber}(:,nn2));
meanArr_perCham = cell(kgroup,1);

for i=1:Kgroup % loop to obtain the average of each cluster per lane 
    meanArr_perCham{i} = nan(Kgroup, length(timeseriesM(1,:)));
    for j=1:CH_no
        ind = ((SPL1{Gnumber}(:,CL_ID) == i) & (chamId == j));
        tmp = timeseriesM(ind,:);
        tmp1 = nanmean(tmp, 1);
        meanArr_perCham{i}(j, :) = tmp1;
    end
end

meanOfChambers = nan(Kgroup, length(meanArr_perCham{2}(1,:)));
for i=1:Kgroup % loop to obtain the mean of each lane
    meanOfChambers(i,:) = nanmean(meanArr_perCham{i}, 1);
end
t=1:(length(meanArr_perCham{2}(1,:)));

%%
f5=figure;

for ik4=group  % loop to plot the average curve for each cluster
    sem = std(meanArr_perCham{ik4}, [], 1, 'omitnan')./sqrt(size(meanArr_perCham{ik4},1));
    hold on
    s1 = shadedErrorBarV2(t, meanOfChambers(ik4,:), 2*sem, 'lineprops', colrs(ik4));
    hold on
    xline(Snumber,'--r') ;%  ,'Starvation'
    xline(Snumber2,'--b') ;%  ,'Starvation'
   ylim([0 limY])
    title(ik4)
    %pause
end
       
   title(['Average' titel1 ]);
   xlim(limX)
   xlabel('Distance from Atg13 foci periphery [pixels]');
   ylabel([titel1 '[a.u.]']); 
   
   
   %    if sae==1
% path_save='C:\Users\aleja\OneDrive\Documents\MATLAB\MATLAB doc\';
% export_path=[path_save titel1 '_ave'];
% hgexport(gcf, export_path, hgexport('factorystyle'), 'Format', 'pdf');
%    else
%    end
   
%     ax = gca;
%     ax.XTickMode = 'manual';
%     xticks(starti:sr:endti);
%     curTick = ax.XTick;
%     ax.XTickLabel = round((curTick)*sr1/60);

% ax = gca;
% ax.XTickMode = 'manual';
% xticks(starti:sr:endti);
% curTick = ax.XTick;
% caxis(axis)
% ax.XTickLabel = round((curTick)*sr1/60);
                        
   
   
   
   if sae==1
saveas(f5,titel1,'pdf')
   else
   end
   
% % % 
% % %  %%%%%%%%----------- obtain heat map ------------------%%%%%%%%%%%%%%%%  
% % %  
% % % groupedMat1=cell(Kgroup,1);
% % % for ik=group % loop to concatenate the clusters on top of each other in a single heatmap
% % % pooledMat_noNan = timeseriesM((SPL1{Gnumber,2}(:,CL_ID) == ik),:);
% % % groupedMat1{ik}= pooledMat_noNan;
% % % end
% % % groupedM1=zeros(size(timeseriesM,1), size(timeseriesM,2));
% % % for ir=1:Kgroup 
% % %     if ir==1
% % %         groupedM1=timeseriesM((SPL1{Gnumber,2}(:,CL_ID) == ir),:);
% % %     else
% % %     groupedM1=[groupedM1;timeseriesM((SPL1{Gnumber,2}(:,CL_ID) == ir),:) ];
% % %     
% % %     end
% % % end
% % % 
% % % 
% % % 
% % % filteredmap = smoothActivityMap(groupedM1, 'SmoothParam', 0.9, 'UpSample', 1);
% % % f6 = figure;imagesc(filteredmap);
% % % 
% % % colorbar;colormap(jet)
% % % xlabel('Time (h)');ylabel('Cell Index')
% % % ax = gca;
% % % ax.XTickMode = 'manual';
% % % xticks(starti:sr:endti);
% % % curTick = ax.XTick;
% % % %caxis(axis)
% % % ax.XTickLabel = round((curTick)*sr1/60);
% % % title(titel1);
% % % 
% % % % % % %    if sae==1
% % % % % % % path_save='C:\Users\aleja\OneDrive\Documents\MATLAB\MATLAB doc\';
% % % % % % % export_path=[path_save titel1 '_heat'];
% % % % % % % hgexport(gcf, export_path, hgexport('factorystyle'), 'Format', 'pdf');
% % % % % % %    else
% % % % % % %    end

           
end 

end

percent=nan(4,2);
for iw=1:2
for ix= 1:4
   
xcells=find(METAD{9}(:,6)==iw & ~isnan(METAD{9}(:,7)))   ;
ycells=find(METAD{9}(:,6)==iw & METAD{9}(:,7)==ix);

percent(ix,iw)=(numel(ycells)./numel(xcells))*100;

end
end

