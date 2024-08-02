
%%AI_RTP: Matlab script to plot the DVH and extract dosimetric parameters
labels_org = {'Spinalcord', 'Esophagus', 'Esophagus_CE','Heart', 'A_Lad','Lung_L','Lung_R','PTV'};
file_path = 'C:\Lulin-home\KBP-lung\CE project\AI_RTP\VCU_Lung_2024_dataset_5';
norg = numel(labels_org);

switch topics

    case 'Plot_DVH'
        %%AI_RTP: Plot DVH for predicted and actual DVHs
        %{
!python train.py --dataroot VCU_Lung_2024_dataset --netG stand_unet --name MAE_Moment_loss --model doseprediction3d --direction AtoB --lambda_L1 1 --dataset_mode dosepred3d --norm batch --batch_size 1 --pool_size 0 --display_id -1 --lr 0.0002 --input_nc 7 --output_nc 1 --display_freq 10 --print_freq 1 --gpu_ids 0
loss1 =  readtable("Z:\LulinY\Lung-dosimetrics\2024\python_code\DoseRTX\results\loss_log.txt", Delimiter={',',':',' ',')'})
w1 = loss1.Var3; w2 = loss1.Var7; w3 = loss1.Var19;

i_tot = w1*58+w2;

figure (1);
plot(i_tot, w3)

   oars = ['lungs','lung_l', 'lung_r','spinalcord', 'esophagus', 'esophagus_ce','heart', 'a_lad',  'ptv']

    oars2 = ['spinalcord', 'esophagus', 'heart',  'ptv']    # for DL
    labels2 = dict.fromkeys(oars,0)
    label_seq = [0, 6, 7, 1, 2, 3, 4, 5, 1]
        %}

        if isfolder(file_path)
            filelist = dir(fullfile(file_path, '*.mat'));  %get list of files and folders in any subfolder
        end
        fl1 = filelist(~[filelist.isdir]);

        color_OAR1={'g-','b-', 'c-','m-','y-','c-','k-','r-'};
        color_OAR2={'g--','b--', 'c--','m--','y--','c--','k--','r--'};
        plotN = 0; fN=0;
        num_case = 60; icase = 0;
        xvar_real = nan(num_case,8,5); xvar_pred= nan(num_case,8,5);
        for i2 = 1: numel(fl1)
            icase = icase+1;
            w5 = fl1(i2); disp(w5.name);
            load(fullfile(file_path, w5.name));


            %% plot DVHs
            plotN = plotN+1;

            index = rem(plotN-1,6)+1;
            %%                plotN = rem(index-1,4)+1;
            if index == 1
                %% figureName = [output_path '/figures/' organ_names_title{igr_org} '-xvar-fit-' num2str(fN) '.jpg'];
                %%print(gcf,'-djpeg', figureName);
                fN = fN+1;
                figure(fN);
                tiledlayout(3,2); nexttile;
                axis([0 70 0 100]); hold on;box on;
            else
                nexttile
            end
            ifh=zeros(norg);
            for iorg = 1:8
                plot(BINS, ...
                    smooth_dvh(HIST_REAL(:,iorg),BINS) *100 , color_OAR1{iorg}, 'LineWidth', 4.); hold on;

                if exist("HIST_PRED")
                    plot(BINS, ...
                        smooth_dvh(HIST_PRED(:,iorg),BINS) *100 , color_OAR2{iorg}, 'LineWidth', 4.);
                end
                xlabel('Dose (Gy)', 'FontSize',14); ylabel ('Volume (%)', 'FontSize',14);
                set(gca,'LineWidth',2,'FontSize',14); grid on;
                legend()
            end



            %{

            for iorg = 1:7

                dose_pres = 1; xaxis_dvh = BINS;
               % disp('org='),iorg
                switch iorg

                    case {6,7}
                        %% for r l lung
                        dv_dose = [5 20];

                        for ic = 1:2
                            d1 = max(2.0,dv_dose(ic))/dose_pres;
                            v1 = interp1(xaxis_dvh,HIST_REAL(:,iorg),d1);
                            xvar_real(icase,iorg, ic) = v1*100;
                        end      %% ic = 1:3

                        ddvh = [diff(HIST_REAL(:,iorg)') 0];
                        dmean = sum(abs(ddvh.*xaxis_dvh));
                        xvar_real(icase,iorg, ic+1) = dmean;
                        ic_plot = [1:3];


                    case {2,4}     %% for esophagus, heart
                        ic = 0;
                        ddvh = [diff(HIST_REAL(:,iorg)') 0];
                        dmean = sum(abs(ddvh.*xaxis_dvh));
                        xvar_real(icase,iorg,ic+1) = dmean;
                        ic_plot = [1];

                    case {3}
                        %% CE
                        dv_dose = [45];

                        for ic = 1:1
                            d1 = max(2.0,dv_dose(ic))/dose_pres;
                            v1 = interp1(xaxis_dvh,HIST_REAL(:,iorg),d1);
                            xvar_real(icase,iorg, ic) = v1*100;
                        end      %% ic = 1:3

                    case {5}
                        %% LAD
                        dv_dose = [15];

                        for ic = 1:1
                            d1 = max(2.0,dv_dose(ic))/dose_pres;
                            v1 = interp1(xaxis_dvh,HIST_REAL(:,iorg),d1);
                            xvar_real(icase,iorg, ic) = v1*100;
                        end      %% ic = 1:3
  case {1}
                                %% spinal

                                for ic = 1:1
                                      v1 = xaxis_dvh(max(find(HIST_REAL(:,iorg)>0.01)));
                                    xvar_pred(icase,iorg, ic) = v1;
                                end      %% ic = 1:3

                end
                        %% for fitted values---------------------------
                        switch iorg

                            case {6,7}
                                %% for r l lung
                                dv_dose = [5 20];
                                for ic = 1:2
                                    d1 = max(2.0,dv_dose(ic))/dose_pres;
                                    v1 = interp1(xaxis_dvh,HIST_PRED(:,iorg),d1);
                                    xvar_pred(icase,iorg, ic) = v1*100;
                                end      %% ic = 1:3

                                ddvh = [diff(HIST_PRED(:,iorg)') 0];
                                dmean = sum(abs(ddvh.*xaxis_dvh));
                                xvar_pred(icase,iorg, ic+1) = dmean;
                                ic_plot = [1:3];


                            case {2,4}     %% for esophagus, heart
                                ic = 0;
                                ddvh = [diff(HIST_PRED(:,iorg)') 0];
                                dmean = sum(abs(ddvh.*xaxis_dvh));
                                xvar_pred(icase,iorg,ic+1) = dmean;
                                ic_plot = [1];

                            case {3}
                                %% CE
                                dv_dose = [45];

                                for ic = 1:1
                                    d1 = max(2.0,dv_dose(ic))/dose_pres;
                                    v1 = interp1(xaxis_dvh,HIST_PRED(:,iorg),d1);
                                    xvar_pred(icase,iorg, ic) = v1*100;
                                end      %% ic = 1:3

                            case {5}
                                %% LAD
                                dv_dose = [15];

                                for ic = 1:1
                                    d1 = max(2.0,dv_dose(ic))/dose_pres;
                                    v1 = interp1(xaxis_dvh,HIST_PRED(:,iorg),d1);
                                    xvar_pred(icase,iorg, ic) = v1*100;
                                end      %% ic = 1:3

  case {1}
                                %% spinal

                                for ic = 1:1
                                      v1 = xaxis_dvh(max(find(HIST_PRED(:,iorg)>0.01)));
                                    xvar_pred(icase,iorg, ic) = v1;
                                end      %% ic = 1:3

                        end %% end switch igr_org

                end
            %}
        end


    case 'MCO'

        xvar_real(:,8,2) = 0.5*(xvar_real(:,6,2) +xvar_real(:,7,2));
        xvar_pred(:,8,2) = 0.5*(xvar_pred(:,6,2) +xvar_pred(:,7,2));
        res = xvar_pred-xvar_real;

        err = nan(8,2);
        for iorg = 1:8
            for ic = 1:2
                err(iorg,ic) = mean(abs(res(:,iorg,ic)),'omitnan');
            end
        end

        tab_err = array2table([err(1,1) err(8,2) err(2, 1) err(3,1) err(4,1)],"VariableNames",{'Spinal Cord Dmax', 'Lung V20', 'Esophagus Dmean', 'Contra Esophagus V20','Heart Dmean'});
        figure(1); hold on
        tiledlayout(1,2); nexttile
        scatter3(res(:,8,2),res(:,2,1), res(:,4,1)); axis([-5 5 -5 5 -5 5])
        set(gca,'FontSize',16);  xlabel('Lung V_{20Gy} (%)', 'FontSize',16); ylabel('Esophagus Mean dose (Gy)', 'FontSize',16);
        zlabel('Heart Mean dose (Gy)', 'FontSize',16);
        nexttile
        scatter3(res(:,8,2),res(:,3,1), res(:,4,1));
        set(gca,'FontSize',16);  xlabel('Lung V_{20Gy} (%)', 'FontSize',16); ylabel('Contra Esophagus Mean dose (Gy)', 'FontSize',16);
        zlabel('Heart Mean dose (Gy)', 'FontSize',16);


    case 'MDS'

        w6 = [res(:,8,2),res(:,2,1), res(:,4,1)];
        Derror = pdist(w6(1:57,:),"euclidean"); [Yerr,eig] = cmdscale(Derror);
        format short g
        [eig cumsum(eig/sum(abs(eig)))]

        meanerr = mean(abs(Derror - pdist(Yerr(:,1:2))),'omitnan')


    case 'dosimetrics'

        %%AI_RTP: Extracted dose metrics
        dose_pres = 60;
        switch igr_org

            case {6,7}
                %% for r l lung
                dv_dose = [5 20];
                xaxis_dvh = BINS;
                %%dv_dose = [0.02*dose_pres(ind_train(index)) 20 30]; %% for
                %% chestwall
                %{
  ic = 1;
  fh2 = @(x) interp1(xaxis_dvh,pcadata_dvh_new(ind_train(index),:,igr_org),x);
    start_value = xaxis_dvh(min(find([(pcadata_dvh_new(ind_train(index),:,igr_org)-0.02)<0])));
    if isempty(start_value)
        start_value = 1.2;
    end
   
    xvar_real(ind_train(index),ic) = fzero(@(x) fh2(x)-0.02, start_value)*100;
                %}
                for ic = 1:2
                    d1 = min(2.0,dv_dose(ic))/dose_pres;
                    v1 = interp1(xaxis_dvh,HIST_REAL(:,igr_org),d1);
                    xvar_real(icase,ic) = v1*100;
                end      %% ic = 1:3

                ddvh = [diff(pcadata_dvh_new(ind_train(index),:,igr_org)) 0];
                dmean = sum(abs(ddvh.*xaxis_dvh));
                xvar_real(ind_train(index),ic+1) = dmean;
                ic_plot = [1:3];


            case 4     %% for esophagus
                dv_dose = [20 60];
                dv_dose = [0.02*dose_pres(ind_train(index)) 20];
                ic = 1;
                fh2 = @(x) interp1(xaxis_dvh,pcadata_dvh_new(ind_train(index),:,igr_org),x);
                start_value = xaxis_dvh(min(find([(pcadata_dvh_new(ind_train(index),:,igr_org)-0.02)<0])));
                if isempty(start_value)
                    start_value = 1.2;
                end

                rel_vol = max(0.03/max(volumePortion_nstd(ind_train(index),3,igr_org),1000),0.02);
                xvar_real(ind_train(index),ic) = fzero(@(x) fh2(x)-rel_vol, start_value)*100;

                for ic = 2:2
                    d1 = min(2.0,dv_dose(ic)/dose_pres(ind_train(index)));
                    v1 = interp1(xaxis_dvh,pcadata_dvh_new(ind_train(index),:,igr_org),d1);
                    xvar_real(ind_train(index),ic) = v1*100;
                end      %% ic = 1:3

                ddvh = [diff(pcadata_dvh_new(ind_train(index),:,igr_org)) 0];
                dmean = sum(abs(ddvh.*xaxis_dvh));
                xvar_real(ind_train(index),ic+1) = dmean;
                ic_plot = [1];

            case {3,17}  %% heart
                ic = 1;
                fh2 = @(x) interp1(xaxis_dvh,pcadata_dvh_new(ind_train(index),:,igr_org),x);
                start_value = xaxis_dvh(min(find([(pcadata_dvh_new(ind_train(index),:,igr_org)-0.02)<0])));
                if isempty(start_value)
                    start_value = 1.2;
                end
                rel_vol = max(0.1/max(volumePortion_nstd(ind_train(index),3,igr_org),1000),0.02);
                xvar_real(ind_train(index),ic) = fzero(@(x) fh2(x)-rel_vol, start_value)*100;

                dv_dose = [0.02 45 50];
                %%         dv_dose = [0.02*dose_pres(ind_train(index)) 50];
                for ic=2:3
                    d1 = min(2.0,dv_dose(ic)/dose_pres(ind_train(index)));
                    v1 = interp1(xaxis_dvh,pcadata_dvh_new(ind_train(index),:,igr_org),d1);
                    xvar_real(ind_train(index),ic) = v1*100;
                end
                ddvh = [diff(pcadata_dvh_new(ind_train(index),:,igr_org)) 0];
                dmean = sum(abs(ddvh.*xaxis_dvh));
                xvar_real(ind_train(index),ic+1) = dmean;



                ic_plot = [1];

            case {5, 9}  %% spinal cord
                ic = 1;
                fh2 = @(x) interp1(xaxis_dvh,pcadata_dvh_new(ind_train(index),:,igr_org),x);
                start_value = xaxis_dvh(min(find([(pcadata_dvh_new(ind_train(index),:,igr_org)-0.5)<0])));
                if isempty(start_value)
                    start_value = 1.2;
                end

                xvar_real(ind_train(index),ic) = fzero(@(x) fh2(x)-0.5, start_value)*100;
                ic = 2;

                start_value = xaxis_dvh(min(find([(pcadata_dvh_new(ind_train(index),:,igr_org)-0.02)<0])));
                if isempty(start_value)
                    start_value = 1.2;
                end

                rel_vol = max(0.03/max(volumePortion_nstd(ind_train(index),3,igr_org),1000),0.02);
                xvar_real(ind_train(index),ic) = fzero(@(x) fh2(x)-rel_vol, start_value)*100;


        end %% end switch igr_org



    case 'output_results'
        %%  a_shape_metrics_scaled = a_shape_metrics;
        fn_stats = "DIBH_stats.csv";
        %{
        for i4= 1:3
            a_shape_metrics_scaled(:,3,i4) = a_shape_metrics(:,3,i4);
        end
        %}
        T21 = cell2table(a_case_ID);
        T22 = array2table([a_shape_metrics_scaled(:,1,1) a_shape_metrics_scaled(:,1,2) a_shape_metrics_scaled(:,1,3)...
            a_shape_metrics_scaled(:,2,1) a_shape_metrics_scaled(:,2,2) a_shape_metrics_scaled(:,2,3)...
            a_shape_metrics_scaled(:,3,1) a_shape_metrics_scaled(:,3,2) a_shape_metrics_scaled(:,3,3)]);
        T2 = [T_case(:,"PatientID") T_case(:,"FB") T_case(:,"C") T_case(:,"B_C") T_case(:,"Site") T22];
        names_metrics = {'Vol Lung_L_FB', 'Vol Lung_R_FB', 'MHD_FB(cm)', 'Vol Lung_L_BH1', 'Vol Lung_R_BH1', 'MHD_BH1(cm)'...
            'Vol Lung_L_BH2', 'Vol Lung_R_BH2', 'MHD_BH2(cm)'};
        names_metrics =  names_metrics([1 4 7 2 5 8 3 6 9]);
        w12 = T2.Properties.VariableNames(1:5);
        T2.Properties.VariableNames = [w12 names_metrics];
        writetable(T2,fn_stats);
        m_vollungl = [a_shape_metrics_scaled(:,1,1) a_shape_metrics_scaled(:,1,2) a_shape_metrics_scaled(:,1,3)];
        m_vollungr = [a_shape_metrics_scaled(:,2,1) a_shape_metrics_scaled(:,2,2) a_shape_metrics_scaled(:,2,3)];
        m_mhd = [a_shape_metrics_scaled(:,3,1) a_shape_metrics_scaled(:,3,2) a_shape_metrics_scaled(:,3,3)];%%%*(-0.1);

        figure(1);
        label_box1 = {'Free Breathing','Chest BH','Belly BH'};
        l_metrics = replace(names_metrics,'_','\_');
        tiledlayout(2,2);
        nexttile;
        h1 = boxplot(m_vollungl,'Labels',label_box1);
        title(l_metrics{1}); grid;
        ylabel(l_metrics{1},'FontSize',16);
        nexttile;
        h1 = boxplot(m_vollungr,'Labels',label_box1);
        title(l_metrics{2}); grid;
        ylabel(l_metrics{2},'FontSize',16);            nexttile;
        h1 = boxplot(m_mhd,'Labels',label_box1);
        title(l_metrics{3}); grid;
        ylabel(l_metrics{3},'FontSize',16);

        print(1,'-djpeg',['DIBH-fig-' int2str(31) '.jpg']);

    case 'paired-analysis'

        %{
        id_rtbrst = [3 6 8 18 22 25 26];
        w12 = isnan(a_shape_metrics_scaled(:,3,1))|isnan(a_shape_metrics_scaled(:,3,2))|isnan(a_shape_metrics_scaled(:,3,3));
        id_MHD = find(~(w12)); id_MHD = id_MHD(1:11);

        a_metrics_new = permute(a_shape_metrics_scaled,[1 3 2]);
        a_metrics_pair = a_metrics_new(id_MHD,:,:);
        %}
        id_rtbrst = [3 6 18 26];
        tiledlayout(2,2);
        nexttile;
        h1 = bar(a_metrics_new(id_MHD,1:3,3)); title('Maximum Heart Distance (cm)'); hold on;
        ylim([-0.5 2]);
        set(gca,'FontSize',16); ylabel('MHD (cm)', 'FontSize',16);
        nexttile;
        h2 = bar(a_metrics_new(id_MHD,1:3,1)); title('Left Lung Volume (CM^3)'); hold on;
        ylim([0 3600]);
        set(gca,'FontSize',16); ylabel('Lt Lung Volume (cm^3)', 'FontSize',16);
        nexttile;

        h3 = bar(a_metrics_new(id_rtbrst,1:3,2)); title('Right Lung Volume (CM^3)'); hold on;
        ylim([0 3600]);

        set(gca,'FontSize',16); ylabel('Rt Lung Volume (cm^3)', 'FontSize',16);

        dv31 = a_metrics_new(id_MHD,1,3) - a_metrics_new(id_MHD,2,3);
        dv32 = a_metrics_new(id_MHD,1,3) - a_metrics_new(id_MHD,3,3);
        dv33 = a_metrics_new(id_MHD,2,3) - a_metrics_new(id_MHD,3,3);

        dv11 = a_metrics_new(id_MHD,1,1) - a_metrics_new(id_MHD,2,1) ;
        dv12 = a_metrics_new(id_MHD,1,1)  - a_metrics_new(id_MHD,3,1) ;
        dv13 = a_metrics_new(id_MHD,2,1)  - a_metrics_new(id_MHD,3,1) ;

        dv21 = a_metrics_new(id_rtbrst,1,2) - a_metrics_new(id_rtbrst,2,2) ;
        dv22 = a_metrics_new(id_rtbrst,1,2)  - a_metrics_new(id_rtbrst,3,2) ;
        dv23 = a_metrics_new(id_rtbrst,2,2)  - a_metrics_new(id_rtbrst,3,2) ;

        mean_dv11 = mean(dv11,'omitnan');
        mean_dv12 = mean(dv12,'omitnan');
        mean_dv13 = mean(dv13,'omitnan');
        mean_dv21 = mean(dv21,'omitnan');
        mean_dv22 = mean(dv22,'omitnan');
        mean_dv23 = mean(dv23,'omitnan');
        mean_dv31 = mean(dv31,'omitnan');
        mean_dv32 = mean(dv32,'omitnan');
        mean_dv33 = mean(dv33,'omitnan');

        a_dv = [mean_dv31 mean_dv11 mean_dv21; mean_dv32 mean_dv12 mean_dv22; mean_dv33 mean_dv13 mean_dv23];



        for ind2 =[1 3]
            pv_pair(ind2,1) = signrank(a_metrics_new(id_MHD,1,ind2),a_metrics_new(id_MHD,2,ind2));
            pv_pair(ind2,2) = signrank(a_metrics_new(id_MHD,1,ind2),a_metrics_new(id_MHD,3,ind2));
            pv_pair(ind2,3) = signrank(a_metrics_new(id_MHD,2,ind2),a_metrics_new(id_MHD,3,ind2));
            for k2 = 1:3
                a_mean(ind2,k2) = mean(a_metrics_new(id_MHD,k2,ind2),'omitnan');
            end

        end
        for ind2 =[2]
            pv_pair(ind2,1) = signrank(a_metrics_new(id_rtbrst,1,ind2),a_metrics_new(id_rtbrst,2,ind2));
            pv_pair(ind2,2) = signrank(a_metrics_new(id_rtbrst,1,ind2),a_metrics_new(id_rtbrst,3,ind2));
            pv_pair(ind2,3) = signrank(a_metrics_new(id_rtbrst,2,ind2),a_metrics_new(id_rtbrst,3,ind2));
            for k2 = 1:3
                a_mean(ind2,k2) = mean(a_metrics_new(id_rtbrst,k2,ind2),'omitnan');
            end
        end

        t_dv = array2table([a_dv(:,1) pv_pair(3,:)' a_dv(:,2) pv_pair(1,:)' a_dv(:,3) pv_pair(2,:)'], ...
            "RowNames",{'FB-DIBH1','FB-DIBH2','DIBH1-DIBH2'},"VariableNames",{'MHD (cm)','pv1','Vol_Lung_L (cc)','pv2','Vol_Lung_R (cc)','pv3'});
        writetable(t_dv,'dv_stats.csv','WriteRowNames',true);


    case 'plot_contour'
        plot3(voxM_PTV(:,1),voxM_PTV(:,2),voxM_PTV(:,3),'b-'); hold on;
        plot3(voxM_Heart(:,1),voxM_Heart(:,2),voxM_Heart(:,3),'r-');
        plot3(voxM_PTV_rot(1,:),voxM_PTV_rot(2,:),voxM_PTV_rot(3,:),'c-');

        figure;
        plot3(voxM_PTV_rot(1,:),voxM_PTV_rot(2,:),voxM_PTV_rot(3,:),'c-'); hold on;
        plot3(surfM(1,:),surfM(2,:),surfM(3,:),'b-'); hold on;
        plot3(structV_rot(1,:),structV_rot(2,:),structV_rot(3,:),'r-');

    case 'plot_dose'

        % Create a 3D plot of the dose distribution
        figure;
        tiledlayout('flow');

        for iz = [16 32 48 64 80 96 112]
            nexttile;
            slice = DOSE_REAL(:, :, iz);
            imagesc(slice);
            colorbar;
            title('Radiation Dose Distribution');
            colormap("jet"); axis equal; view(2)
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        end

        figure;
        tiledlayout('flow');

        for iz = [16 32 48 64 80 96 112]
            nexttile;
            slice = DOSE_PRED(:, :, iz);
            imagesc(slice);
            colorbar;
            title('Radiation Dose Distribution');
            colormap("jet"); axis equal; view(2)
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        end

end    %% topics

