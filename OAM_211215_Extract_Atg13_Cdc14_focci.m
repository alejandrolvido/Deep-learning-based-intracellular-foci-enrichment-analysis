


% to extract nuclear translocators in all channels
parpool('local',2)

type='.tif';
path_h='C:\Users\aleja\OneDrive\Documents\MATLAB\MATLAB_doc\';

exp_name='211124_cdc14_atg13_foci'; 
positions=0;%:2

disp(positions)   
% TFP GFP/mNG mKok mRuby mNeptune CYOFP
fl_used=[0 1 0 1 0 0 0];
% =========================================================================
bin_used=1;
% W5_mean_modifier=1;
% IT_mean_modifier= 1;
% G_mean_modifier= 1;
% K_mean_modifier= 1; 
% R_mean_modifier= 1;
% FR_mean_modifier= 1;
% vac_mean_modifier=1;
% E6h_mean_modifier=1;
% peak_cutoff=0.75;

max_size_vs_largest_cell=3;
atg13_mean_modifier=3.0;

tic
pos=0;
for ptmp_h=1:405 % loop to extract everyframe in a movie
    
 
%     path_seg='D:\N_folder\new_segmentations\';
     path_seg='C:\Users\aleja\OneDrive\Documents\MATLAB\MATLAB_doc\';
     Imagepath=[path_seg exp_name '\Pos' num2str(pos) '\'];
    %path_seg='C:\Users\aleja\OneDrive\Documents\MATLAB doc\Scores_Home\';
%     load([path_seg exp_name '\Seeds\img_000000' num2str(ptmp_h) '_BrightField_000_cp_masks'])
    
       image_number=sprintf('%09d',ptmp_h);        
       Lcells=imread([path_seg exp_name '\Seeds\img_' image_number '_BrightField_000_cp_masks.png']);
%     imagesc(Lcells)
      
    numbM=1;
    x_size =size(Lcells,1);   %1040;
    y_size =size(Lcells,2);   %1388;
    
    appr_vol                      =zeros((max(Lcells(:))),numbM);
    %---------------------------------------------------------------------
    %Allocation of features extracted from all fluorophores regardless of the
    %protein tagged

           cell_suma_mNG=cell((max(Lcells(:))),1);
           cell_averag_mNG=cell((max(Lcells(:))),1);
           cell_maxima_mNG=cell((max(Lcells(:))),1);
           cell_variation_mNG=cell((max(Lcells(:))),1);  

           cell_suma_Ru=cell((max(Lcells(:))),1);
           cell_averag_Ru=cell((max(Lcells(:))),1);
           cell_maxima_Ru=cell((max(Lcells(:))),1);
           cell_variation_Ru=cell((max(Lcells(:))),1);    
 
            cell_Total_mNG=cell((max(Lcells(:))),1);  
            cell_Total_Ru=cell((max(Lcells(:))),1);  
           
% %Neon_Green==================================================================================
%     nucl_area_Neon_Green                 =zeros((max(Lcells(:))),numbM);
%     cyt_area_Neon_Green                  =zeros((max(Lcells(:))),numbM);
%     mean_Neon_Green_int_per_area_C        =zeros((max(Lcells(:))),numbM);
%     mean_Neon_Green_int_per_area_N        =zeros((max(Lcells(:))),numbM);
%     Neon_Green_Conc_T                     =zeros((max(Lcells(:))),numbM);
%     Neon_Green_Conc_C                     =zeros((max(Lcells(:))),numbM);
%     Neon_Green_Conc_N                     =zeros((max(Lcells(:))),numbM);
%     nucl_Neon_Green                      =false(x_size,y_size,numbM);
%     Neon_Green_mean_int_N               =zeros((max(Lcells(:))),numbM);
%     Neon_Green_mean_int_C               =zeros((max(Lcells(:))),numbM);
%     Neon_Green_loc                        =false(x_size,y_size,numbM);
%     tot_Neon_Green_fl_cyt = zeros((max(Lcells(:))),numbM);
%     tot_Neon_Green_fl_nuc = zeros((max(Lcells(:))),numbM);
% 
%     %Ruby==================================================================================
%     nucl_area_Ruby                  =zeros((max(Lcells(:))),numbM);
%     cyt_area_Ruby                   =zeros((max(Lcells(:))),numbM);
%     mean_Ruby_int_per_area_C        =zeros((max(Lcells(:))),numbM);
%     mean_Ruby_int_per_area_N        =zeros((max(Lcells(:))),numbM);
%     Ruby_Conc_T                     =zeros((max(Lcells(:))),numbM);
%     Ruby_Conc_C                     =zeros((max(Lcells(:))),numbM);
%     Ruby_Conc_N                     =zeros((max(Lcells(:))),numbM);
%     nucl_Ruby                       =false(x_size,y_size,numbM);
%     Ruby_mean_int_N                 =zeros((max(Lcells(:))),numbM);
%     Ruby_mean_int_C                 =zeros((max(Lcells(:))),numbM);
%     Ruby_loc                        =false(x_size,y_size,numbM);
%     tot_Ruby_fl_cyt                 = zeros((max(Lcells(:))),numbM);
% %     tot_Ruby_fl_nuc                 = zeros((max(Lcells(:))),numbM);
%     
%  mean_T_int_per_area_Ru  =zeros((max(Lcells(:))),numbM);
%  T_Conc_Ru               =zeros((max(Lcells(:))),numbM);
    
    tmp_sizes=zeros(1,(max(Lcells(:)))); %allocation
    Lcells1=Lcells;
    %par
    for i5=1:(max(Lcells(:)))
        tmp_sizes(i5)=tmp_sizes(i5)+sum(sum(Lcells1==i5));
    end
    max_allowed_cell_size=max_size_vs_largest_cell*max(tmp_sizes);
    s=round(sqrt(max_allowed_cell_size))+50; %this will determine the size of the matrix we will put put_vac etc into
    
%     for c_time=1;ptmp_h%;1:(max(Lcells(:)))
%         current_time=c_time;  
%         disp(c_time)
%         Lcells=all_obj.cells(:,:,c_time); %a broadcast variable
        
        %Variables to be used in construction of V_loc E6_loc nucl_mC W5_loc nucl_K
        V_loc_tmp=zeros(s,s,(max(Lcells(:))));
        location_cell=zeros((max(Lcells(:))),4);
      
%         nucl_IT_tmp=zeros(x_size,y_size);
%         nucl_K_tmp=zeros(x_size,y_size);
        nucl_Green_tmp=zeros(x_size,y_size);
%         nucl_Nep_tmp=zeros(x_size,y_size);
        nucl_Ru_tmp=zeros(x_size,y_size);
        
        w5_loc_tmp=zeros(s,s,(max(Lcells(:))));
        nuc_mC_tmp=zeros(x_size,y_size);
        %Find Median Fluorescence
        image_number=sprintf('%09d',ptmp_h);

        if fl_used(1,2)==1
%             % IT=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_470 mTFP_000.tif']);
            IG=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_505 mNG_000.tif']);
%              IK=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_555 new_mKok_000.tif']);
%            % IK=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_555 mKok_000.tif']);
%             ICa=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_615 nm_000.tif']);
            
            IG=double(IG);
            IG=medfilt2(IG,'symmetric');
            backgr_G=(((IG)+1).*(~Lcells));
            backgr_G=sort(backgr_G(:));
            [tmpV,posH]=max(backgr_G>0);
            backgr_G=median(backgr_G(max([1 posH-1]):end))-1;
            IG=IG-backgr_G;
            all_obj.med_G(1,ptmp_h)  =backgr_G;
        end

        if fl_used(1,4)==1
           % IRu=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_555 new_mRuby3_000.tif']);
           IRu=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_555 mRuby3_000.tif']);
            IRu=double(IRu);
            IRu=medfilt2(IRu,'symmetric');
            backgr_Ru=((double(IRu)+1).*(~Lcells));
            backgr_Ru=sort(backgr_Ru(:));
            [~,posH]=max(backgr_Ru>0);
            backgr_Ru=median(backgr_Ru(max([1 posH-1]):end))-1;
            IRu=(IRu-backgr_Ru);
            all_obj.med_Ru(1,ptmp_h)  =backgr_Ru;
        end

        %par
        for cell_no=1:(max(Lcells(:)))
           % if c_time==ptmp_h%>cell_exists(cell_no,2)
                
                ccell=(Lcells==cell_no); %a temporary variable imagesc(ccell)
                appr_vol(cell_no)  = appr_vol(cell_no)+(sum(ccell(:))).^(1.5); % 2D area ^(3/2) %reduction variable
                
                %Get Cell Location imagesc(ccell) imagesc(Lcells)
                cell_margin=1;
                [x_cn,y_cn]=get_wind_coord1_EZGI(ccell,cell_margin);
                location_cell(cell_no,:)=location_cell(cell_no,:)+[y_cn(1) y_cn(end) x_cn(1) x_cn(end)];             

                % define the cell in each channel
                   cell0=double(ccell(y_cn,x_cn)); % whole cell
                   %imagesc(cell0)
                   cell00=bwmorph(cell0,'remove'); %periphery cell
                   %imagesc(cell00)
                    ccell2=double(ccell(y_cn,x_cn)).*double(IRu(y_cn,x_cn)); % imagesc(ccell2) % cell in ruby3
                    ccell3=double(ccell(y_cn,x_cn)).*double(IG(y_cn,x_cn));  % imagesc(ccell3) % cell in mNG
                   
                    cell_Total_Ru{cell_no}=sum(ccell2(:));
                    cell_Total_mNG{cell_no}=sum(ccell3(:));
                    
 % define foci based on mean fluorescence intensity        atg13_mean_modifier=2           
                    put_atg12=(ccell2>=(atg13_mean_modifier.*mean(ccell2(ccell2~=0))));
                   % figure;imagesc(put_atg12)
                   
                   put_atg13=bwlabel(put_atg12);
                   % imagesc(put_atg13)

                 if isempty(put_atg13)==0 && sum(put_atg13(:))~=0
                     
                      suma_Ru=nan(10,20);
                      averag_Ru=nan(10,20);
                      maxima_Ru=nan(10,20);
                      variation_Ru=nan(10,20);  
                      
                      suma_mNG=nan(10,20);
                      averag_mNG=nan(10,20);
                      maxima_mNG=nan(10,20);
                      variation_mNG=nan(10,20);    
                      
                   for iyt=1:max(put_atg13(:)) 
                     
                        put_atg=(put_atg13==iyt) ; % imagesc(put_atg)  
                        stopfact=0;
                        A=1;
%                         for iiu=1:20
                       while stopfact~=2
                     
                      put_at=imdilate(put_atg,[1 1 1; 1 1 1; 1 1 1]);%
                      put_at1=(put_at-put_atg);
                      put_at2=ccell2.*put_at1;
                      put_at3=ccell3.*put_at1;
                         
                      suma_mNG(iyt,A)=sum(put_at3(put_at3(:)~=0));
                      averag_mNG(iyt,A)=nanmean(put_at3(put_at3(:)~=0)); % imagesc(put_at3)
                      
                     if ~isempty(max(put_at3(put_at3(:)==0)))
                     maxima_mNG(iyt,A)=nan;  
                     else
                       maxima_mNG(iyt,A)=max(put_at3(put_at3(:)~=0)); 
                     end
                     
                      variation_mNG(iyt,A)=nanstd(put_at3(put_at3(:)~=0)); 
                      
                      suma_Ru(iyt,A)=sum(put_at2(put_at2(:)~=0));
                      averag_Ru(iyt,A)=nanmean(put_at2(put_at2(:)~=0));
                       
                      if ~isempty(max(put_at2(put_at2(:)==0)))
                      maxima_Ru(iyt,A)=nan; 
                     else
                      maxima_Ru(iyt,A)=max(put_at2(put_at2(:)~=0));  
                     end
                      
%                       maxima_Ru(iyt,A)=max(put_at2(put_at2(:)~=0));
                      variation_Ru(iyt,A)=nanstd(put_at2(put_at2(:)~=0));
                      
%                        figure; imagesc(put_at2)
%                        figure; imagesc(put_at3)
                       put_atcell=put_at+cell00;
%                        imagesc(put_atcell)
                       stopfact=max(put_atcell(:));
                       put_atg=put_at;
%                        pause 
%                        close all
                       A=A+1;
                       end

                   end

                      cell_suma_mNG{cell_no}=suma_mNG;
                      cell_averag_mNG{cell_no}=averag_mNG;
                      cell_maxima_mNG{cell_no}=maxima_mNG;
                      cell_variation_mNG{cell_no}=variation_mNG;   
     
                      cell_suma_Ru{cell_no}=suma_Ru;
                      cell_averag_Ru{cell_no}=averag_Ru;
                      cell_maxima_Ru{cell_no}=maxima_Ru;
                      cell_variation_Ru{cell_no}=variation_Ru;  
                      
                 end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mTFP MTFP
    
               
           % end %end if cell exists
        end %parfor
%     end %time-loop
    
    %add all the fields to all_obj structure
% %     all_obj.appr_vol =appr_vol;
% %     
% %     %------------------------------------TFP----------------------------
% %     if fl_used(1,1)==1
% %         all_obj.mean_T_int_per_area_T        =mean_IT_int_per_area_T;
% %         all_obj.T_Conc_T                     =IT_Conc_T;
% %         all_obj.max_nucl_int_IT              =max_nucl_int_IT;
% %            
% %     end
% %     %------------  NG ------------------------------------------------
% %     if fl_used(1,2)==1
% %         all_obj.mean_Green_int_per_area_T    =mean_Green_int_per_area_T;
% %         all_obj.Green_Conc_T                 =Green_Conc_T;
% %          all_obj.max_nucl_int_Green          =max_nucl_int_Green;   
% %     end
% %     %------------------ mKok ---------------------------------------------------
% %     if fl_used(1,3)==1
% %         all_obj.mean_K_int_per_area_T        =mean_K_int_per_area_T;
% %         all_obj.K_Conc_T                     =K_Conc_T;
% %         all_obj.max_nucl_int_K               =max_nucl_int_K;
% %     end
% % 
% %     %--------------------------555 Ruby------------------------------------------
% %     if fl_used(1,4)==1
% %                       
% %         all_obj.mean_Ru_int_per_area_T      =mean_Ru_int_per_area_T;        
% %         all_obj.Ru_Conc_T                   =Ru_Conc_T;                   
% %         all_obj.max_nucl_int_Ru             =max_nucl_int_Ru;          
% %         
% %     end
% %     %---------------------------------------------------------------------
% % 
% %     %--------------------------615 nm------------------------------------------
% %     if fl_used(1,5)==1
% %                        
% %         all_obj.mean_Nep_int_per_area_T     =mean_Nep_int_per_area_T;
% %         all_obj.Nept_Conc_T                 =Nept_Conc_T;
% %         all_obj.max_nucl_int_Nep            =max_nucl_int_Nep;
% %         
% %     end
% %    %-----------------------------------------------------------------------------
% %    
% %     if fl_used(1,6)==1
% %          all_obj.Cdc10_STD     =Cdc10_stds;
% %          all_obj.Cdc10_Means   =Cdc10_means;
% % %         all_obj.Cdc10_status  =Cdc10_status;
% % 
% %     end

%     %%%Neon_Green
%     all_obj.nucl_area_Neon_Green                  =nucl_area_Neon_Green;
%     all_obj.cyt_area_Neon_Green                   =cyt_area_Neon_Green;
%     all_obj.mean_Neon_Green_int_per_area_C        =mean_Neon_Green_int_per_area_C ;
%     all_obj.mean_Neon_Green_int_per_area_N        =mean_Neon_Green_int_per_area_N;
%     all_obj.Neon_Green_Conc_T                     =Neon_Green_Conc_T;
%     all_obj.Neon_Green_Conc_C                     =Neon_Green_Conc_C;
%     all_obj.Neon_Green_Conc_N                     =Neon_Green_Conc_N;
%     all_obj.nucl_Neon_Green                       =nucl_Neon_Green ;
%     all_obj.Neon_Green_mean_int_N                 =Neon_Green_mean_int_N   ;
%     all_obj.Neon_Green_mean_int_C                 =Neon_Green_mean_int_C ;
%     all_obj.Neon_Green_loc                        =Neon_Green_loc ;
%     all_obj.tot_Neon_Green_fl_cyt                 =tot_Neon_Green_fl_cyt;
%     all_obj.tot_Neon_Green_fl_nuc                 =tot_Neon_Green_fl_nuc;
%     
%     %%%Ruby
%     all_obj.nucl_area_Ruby                  =nucl_area_Ruby;
%     all_obj.cyt_area_Ruby                   =cyt_area_Ruby;
%     all_obj.mean_Ruby_int_per_area_C        =mean_Ruby_int_per_area_C ;
%     all_obj.mean_Ruby_int_per_area_N        =mean_Ruby_int_per_area_N;
%     all_obj.Ruby_Conc_T                     =Ruby_Conc_T;
%     all_obj.Ruby_Conc_C                     =Ruby_Conc_C;
%     all_obj.Ruby_Conc_N                     =Ruby_Conc_N;
%     all_obj.nucl_Ruby                       =nucl_Ruby ;
%     all_obj.Ruby_mean_int_N                 =Ruby_mean_int_N   ;
%     all_obj.Ruby_mean_int_C                 =Ruby_mean_int_C ;
%     all_obj.Ruby_loc                        =Ruby_loc ;
%     all_obj.tot_Ruby_fl_cyt                 =tot_Ruby_fl_cyt;
%     all_obj.tot_Ruby_fl_nuc                 =tot_Ruby_fl_nuc;
           
%      'cell_suma_mNG','cell_averag_mNG', 'cell_maxima_mNG', 'cell_variation_mNG','cell_suma_Ru','cell_averag_Ru', 'cell_maxima_Ru', 'cell_variation_Ru',    
                                            
   path_save='C:\Users\aleja\OneDrive\Documents\MATLAB\MATLAB_doc\Scores_Home\';
   name1=[exp_name '_pos' num2str(ptmp_h) '_extr_gradient'];
   % name1=[exp_name '/' exp_name '_pos' num2str(pos) '_extr_fl_data_Cdc10'];
    save(fullfile(path_save,name1),'atg13_mean_modifier','cell_suma_mNG','cell_averag_mNG', 'cell_maxima_mNG', 'cell_variation_mNG','cell_suma_Ru','cell_averag_Ru', 'cell_maxima_Ru', 'cell_variation_Ru','cell_Total_mNG','cell_Total_Ru');

      
                    
    
end
toc
 delete(gcp('nocreate'))


