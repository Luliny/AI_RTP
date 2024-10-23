
%%AI_RTP: Matlab script to plot the DVH and extract dosimetric parameters
labels_org = {'Spinalcord', 'Esophagus', 'Esophagus_CE','Heart', 'A_Lad','Lung_L','Lung_R','PTV'};
file_path_info = 'C:\Lulin-home\KBP-lung\CE project\AI_RTP\VCU_Lung_2024_dataset_5';
file_path_info = 'Z:\LulinY\Lung-dosimetrics\2024\python_code\AI_RTP';
file_path = 'R:\LulinY\Processed_Lung_Data_2';
file_path = 'R:\LulinY\Processed_NPZ';
figure_path = 'Z:\LulinY\Lung-dosimetrics\2024\python_code\AI_RTP\figures';
pat_case_id = "VCU_Lung_"+digitsPattern(3);
norg = numel(labels_org);

switch topics

    case 'plot_DVH'
        %%AI_RTP: Plot DVH for predicted and actual DVHs
     

        if isfolder(file_path)
            filelist = dir(fullfile(file_path, '*.mat'));  %get list of files and folders in any subfolder
        end
        fl1 = filelist(~[filelist.isdir]);

        %% read case info table
        file_case_info = fullfile(file_path_info, 'case_info_structures_25June2024.csv');
        case_info = readtable(file_case_info,"ReadRowNames",true);

        
        color_OAR1={'g-','b-', 'c-','m-','y-','k--','k-','r-'};
        color_OAR2={'g--','b--', 'c--','m--','y--','k--','k--','r--'};
        plotN = 0; fN=0;
        num_case = 100; icase = 0;

        for i2 = 1:numel(fl1)
            icase = icase+1;
            w5 = fl1(i2); disp(w5.name);
            load(fullfile(file_path, w5.name));
            w6 = w5.name; [dummy, w7, ext] = fileparts(w6);  %%strlength(w6)-9);
            id_case_tmp = extract(w7,pat_case_id);
            id_case1 = id_case_tmp{1};  id_case1 = convertCharsToStrings(id_case1); 


            dose_pres = case_info{id_case1,"PrescripedDose_cGy_"}/100; xaxis_dvh = BINS/dose_pres;


            %% plot DVHs
            plotN = plotN+1;

            index = rem(plotN-1,9)+1;
            %%                plotN = rem(index-1,4)+1;
            if index == 1

                fN = fN+1;
                figure(fN);
                tiledlayout(3,3,'TileSpacing','compact'); nexttile;
                axis([0 70 0 100]); hold on;box on; title(strrep(w7,'_','-'),'FontSize',6)
            else
                nexttile; axis([0 70 0 100]); hold on;box on; title(strrep(w7,'_','-'), 'FontSize',6)
            end
            ifh=zeros(1,norg);
            for iorg = 1:8
               ifh(iorg)= plot(BINS, ...
                    smooth_dvh(HIST_REAL(:,iorg),BINS) *100 , color_OAR1{iorg}, 'LineWidth', 2.); hold on;

                if exist("HIST_PRED")
                    plot(BINS, ...
                        smooth_dvh(HIST_PRED(:,iorg),BINS) *100 , color_OAR2{iorg}, 'LineWidth', 2.);
                end
                xlabel('Dose (Gy)', 'FontSize',6); ylabel ('Volume (%)', 'FontSize',6);
                set(gca,'LineWidth',2,'FontSize',6); grid on;               
            end

           if index == 1 && plotN>2
            lgd = legend(ifh,labels_org); lgd.Layout.Tile = 1;
            if fN>1.5
               %% figureName = [figure_path '/figures/' organ_names_title{igr_org} '-xvar-fit-' num2str(fN) '.jpg'];
               figureName = fullfile(figure_path, ['fig-' num2str(fN) '-DVH.jpg']);
                print(fN-1,'-djpeg', figureName);
   %% uncomment this to save figures
            end
           end
        end      %% cases


    case 'dosimetrics'

        if isfolder(file_path)
            filelist = dir(fullfile(file_path, '*.mat'));  %get list of files and folders in any subfolder
        end
        fl1 = filelist(~[filelist.isdir]);

        %% read case info table
        file_case_info = fullfile(file_path_info, 'case_info_structures_25June2024.csv');
        case_info = readtable(file_case_info,"ReadRowNames",true);
        file_plan_info = 'R:\LulinY\Processed_NPZ\Cases_processing_summary_new.csv';
        plan_info = readtable(file_plan_info,"ReadRowNames",true);        
        num_pl = size(plan_info,1); iplan  = 0;
%%        Dx = nan(num_pl,1);
        for i2 = 1: num_pl
            iplan = iplan+1;
            id_case = plan_info.CaseID{i2};
            id_case_tmp = extract(id_case,pat_case_id);
            id_case1 = id_case_tmp{1};  id_case1 = convertCharsToStrings(id_case1); 
            dose_pres = case_info{id_case1,"PrescripedDose_cGy_"}/100;
            plan_info.Dx(i2) = dose_pres;
        end
        writetable(plan_info,'plan_info_dx.csv');

        return

        color_OAR1={'g-','b-', 'c-','m-','y-','c-','k-','r-'};
        color_OAR2={'g--','b--', 'c--','m--','y--','c--','k--','r--'};
        plotN = 0; fN=0;
        num_case = 100; icase = 0;
        xvar_real = nan(num_case,9,10); %%xvar_pred= nan(num_case,9,10);
        for i2 = 1: 30  %%%numel(fl1)
            icase = icase+1;
            w5 = fl1(i2); disp(w5.name);
            load(fullfile(file_path, w5.name));
            w6 = w5.name; [dummy, w7, ext] = fileparts(w6);  %%strlength(w6)-9);
            id_case_tmp = extract(w7,pat_case_id);
            id_case1 = id_case_tmp{1};  id_case1 = convertCharsToStrings(id_case1); 


            dose_pres = case_info{id_case1,"PrescripedDose_cGy_"}/100; xaxis_dvh = BINS/dose_pres;


     %%  labels_org = {'Spinalcord', 'Esophagus', 'Esophagus_CE','Heart', 'A_Lad','Lung_L','Lung_R','PTV'};    
            for iorg = 1:7

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
                        xvar_real(icase,iorg, ic+1) = dmean*dose_pres;
                        ic_plot = [1:3];


                    case {2}     %% for esophagus

                            ic = 1;
            fh2 = @(x) interp1(xaxis_dvh,smooth_dvh(HIST_REAL(:,iorg),xaxis_dvh),x);
            start_value = xaxis_dvh(min(find((smooth_dvh(HIST_REAL(:,iorg),xaxis_dvh)-0.5)<0)));
    if isempty(start_value)
        start_value = 1.2;
    end
xvar_real(icase,iorg,ic) = fzero(@(x) fh2(x)-0.01, start_value)*dose_pres;

                        dv_dose = [0 30 60];

                        for ic = 2:3
                            d1 = max(2.0,dv_dose(ic))/dose_pres;
                            v1 = interp1(xaxis_dvh,HIST_REAL(:,iorg),d1);
                            xvar_real(icase,iorg, ic) = v1*100;
                        end      %% ic = 1:3

                        ic = 4;
                        ddvh = [diff(HIST_REAL(:,iorg)') 0];
                        dmean = sum(abs(ddvh.*xaxis_dvh));
                        xvar_real(icase,iorg,ic+1) = dmean*dose_pres;
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
%{
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
            end   %% for each organ
        end    %% for each case

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




    case "plot_contours_nrrd"
        CaseID = 'VCU_Lung_006';

        if isfolder(file_path)
            filelist1 = dir(fullfile(file_path, [CaseID '*CT.nrrd']));  %get list of files and folders in any subfolder
        filelist2 = dir(fullfile(file_path, [CaseID '*PTV.nrrd']));
        filelist3 = dir(fullfile(file_path, [CaseID '*RTSTRUCTS.nrrd']));
        filelist4 = dir(fullfile(file_path, [CaseID '*dose.nrrd']));
       % filelist5 = dir(fullfile(file_path, [CaseID '*.mat']));
        end

        CT1 = nrrdread(fullfile(file_path, filelist1.name));
        PTV = nrrdread(fullfile(file_path, filelist2.name));
        OAR = nrrdread(fullfile(file_path, filelist3.name));
        DOSE = nrrdread(fullfile(file_path, filelist4.name));
       % PlanData = load(fullfile(file_path, filelist5.name));

      % CT_mat = PlanData.CT;
%%vol1 = volshow(CT1);
vol1 = sliceViewer(CT1);
vol2 = volshow(OAR);
vol4 = volshow(PTV);
vol3 = sliceViewer(DOSE,"DisplayRangeInteraction","on", "Colormap",turbo);
colorbar;

%%vol6 = volshow(CT_mat);

    case "plot_contours_mat"
        CaseID = 'VCU_Lung_103_0';

       if isfolder(file_path)
        filelist5 = dir(fullfile(file_path, [CaseID '*.mat']));
       end

      PlanData = load(fullfile(file_path, filelist5.name));

       CT_mat = PlanData.CT;
       dose_mat = PlanData.DOSE_REAL;
       oar_mat = PlanData.OAR;
       ptv_mat = PlanData.PTV;
%%vol1 = volshow(CT1);
figure(1);
vol1 = sliceViewer(CT_mat);
figure(2);
vol3 = sliceViewer(dose_mat,"DisplayRangeInteraction","on", "Colormap",turbo);
colorbar;
figure(3);
vol4 = volshow(oar_mat);
figure(4);
vol5 = volshow(ptv_mat);

figure(5);
vol4 = sliceViewer(oar_mat,"DisplayRangeInteraction","on", "Colormap",turbo);
%%vol6 = volshow(CT_mat);


    case "plot_contours_mat_reduced"
        CaseID = 'VCU_Lung_103_0.npz_reduced';

       if isfolder(file_path)
        filelist5 = dir(fullfile(file_path, [CaseID '*.mat']));
       end

      PlanData = load(fullfile(file_path, filelist5.name));

       %CT_mat = PlanData.CT;
       dose_mat = PlanData.DOSE;
       oar_mat = PlanData.OAR;
       %ptv_mat = PlanData.PTV;
%%vol1 = volshow(CT1);
%figure(1);
%vol1 = sliceViewer(CT_mat);
figure(2);
vol3 = sliceViewer(dose_mat,"DisplayRangeInteraction","on", "Colormap",turbo);
colorbar;
figure(3);
vol4 = volshow(oar_mat);
%figure(4);
%vol5 = volshow(ptv_mat);

figure(5);
vol4 = sliceViewer(oar_mat,"DisplayRangeInteraction","on", "Colormap",turbo);
%%vol6 = volshow(CT_mat);

end    %% topics

