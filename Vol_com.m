clc, close all, clear all
%%
figure
subplot(1,2,1)
[x,yp] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P1_Intra-OP - Preop- DVH Data.txt");

x(1,1) = 0;
yp(1,1) = 0; 

plot(x,yp,'r--')
hold on
[x,yp_io] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P1_Intra-OP - Intraop - DVH Data.txt");

x(1,1) = 0;
yp_io(1,1) = 0;
plot(x,yp_io,'r')

% OAR 1
[x,yr] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P1_Rectum Intra-OP - Preplan - DVH Data.txt");
plot(x,yr,'g--')

[x,yr_io] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P1_Rectum Intra-OP - Intraop - DVH Data.txt");
plot(x,yr_io,'g')

% OAR 2
[x,yu] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P1_Urethra Intra-OP - Preplan - DVH Data.txt");
plot(x,yu,'b--')

[x,yu_io] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P1_Urethra Intra-OP - Intraop - DVH Data.txt");
plot(x,yu_io,'b')

legend("Plan Target","IntraOp Target","Plan Rectum","IntraOp Rectum","Plan Urethra","IntraOp Urethra")
xlabel('Dose (Gy)')
ylabel('% Volume')

subplot(1,2,2)
[x,yp_p2] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P2_Intra-OP - Preop- DVH Data.txt");

x(1,1) = 0;
yp_p2(1,1) = 0; 

plot(x,yp_p2,'r--')
hold on

[x,yp_io_p2] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P2_Intra-OP - Intraop - DVH Data.txt");
x(1,1) = 0;
yp_io_p2(1,1) = 0; 

plot(x,yp_io_p2,'r')
legend("Plan Target_ P2","IntraOp Target_ P2")
xlabel('Dose (Gy)')
ylabel('% Volume')
%% Patient 2

figure
subplot(1,2,1)
[x,yp_p2] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P2_Intra-OP - Preop- DVH Data.txt");

x(1,1) = 0;
yp_p2(1,1) = 0; 

plot(x,yp_p2,'r--')
hold on
[x,yp_io_p2] = PullDVH("C:\Users\SA\Desktop\0531\Thesis\Data\histogram_DVH_pre and Intra op plans\P2_Intra-OP - Intraop - DVH Data.txt");


%% Find Volume difference
% https://www.mathworks.com/help/matlab/math/integration-of-numeric-data.html

%% Delta_vol_p1
xverts_plan_p1 = x(1:1:end);
yverts_plan_p1 = yp(1:1:end);
p1 = patch(xverts_plan_p1,yverts_plan_p1,'g');%,'LineWidth',0,5);
Vol_plan_p1 = trapz(yp);
cvol_plan_p1 = cumtrapz(yp);
% T = table(x',cvol_plan','VariableNames',{'Dose','CumulativeVol'})

yverts_io = [yp_io(1:end)];
Vol_io_p1 = trapz(yp_io);
cvol_io_p1 = cumtrapz(yp_io);
Delta_vol_p1 = Vol_io_p1 - Vol_plan_p1

%% Delta_vol_p2
xverts_plan_p2 = x(1:1:end);
yverts_plan_p2 = yp_p2(1:1:end);
p2 = patch(xverts_plan_p2,yverts_plan_p2,'w');%,'LineWidth',0,5);
Vol_plan_p2 = trapz(yp_p2);
cvol_plan_p2 = cumtrapz(yp_p2);
% T = table(x',cvol_plan','VariableNames',{'Dose','CumulativeVol'})

yverts_io_p2 = [yp_io_p2(1:end)];
Vol_io_p2 = trapz(yp_io_p2);
cvol_io_p2 = cumtrapz(yp_io_p2);
Delta_vol_p2 = Vol_io_p2 - Vol_plan_p2

subplot(1,2,2)
plot(cvol_plan_p1,'r')
hold on
plot(cvol_io_p1,'r--')
hold on
plot(cvol_plan_p2,'b')
hold on
plot(cvol_io_p2,'b--')

legend("Prostate Plan Vol_ P1","Prostate IntraOp Vol_ P1","Prostate Plan Vol_ P2","Prostate IntraOp Vol_ P2",'Location','southeast')
title('Cumulative vol')
xlabel('Dose (Gy)')
ylabel('Vol')
%% P1 pre
dose = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_Preplan\DO001.dcm'));
doseI = dicominfo('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_Preplan\DO001.dcm');
total_area_pre = 0;
col = [];
for i = 1:size(dose,3)
    figure;
    imshow(dose(:,:,i),[]);axis on; hold on
    contour(dose(:,:,i)*doseI.DoseGridScaling,'ShowText','on');set(gca, 'YDir','reverse');
    if sum((dose(:,:,i)*doseI.DoseGridScaling)>=190,'all') >=1
        hold on
        contour(dose(:,:,i)*doseI.DoseGridScaling,[190,190],'LineColor','r','LineWidth',3,'ShowText','on');set(gca, 'YDir','reverse');
        [area,~,~] = Contour2Area(contour(dose(:,:,i)*doseI.DoseGridScaling,[190,190]));
        title(['area: ',num2str(area)])
%         total_area_pre = total_area_pre + area;
        col = [col;[i, sum(area,'all')]];
    end
%     waitforbuttonpress
    close
end 

for i = 1:length(col)-1
    volu_inter = (col(i,2) + col(i+1,2))*5/2;
    total_area_pre = total_area_pre + volu_inter;
end 
disp(sum(total_area_pre/1000,'all'))

%% P1 interop
dose = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_IntraoperativeTrementPlan\DO001.dcm'));
doseI = dicominfo('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_IntraoperativeTrementPlan\DO001.dcm');
total_area_inter = 0;
col = [];
for i = 1:size(dose,3)
    figure;
    imshow(dose(:,:,i),[]);axis on; hold on
    contour(dose(:,:,i)*doseI.DoseGridScaling,'ShowText','on');set(gca, 'YDir','reverse');
    if sum((dose(:,:,i)*doseI.DoseGridScaling)>=180,'all') >=1
        hold on
        contour(dose(:,:,i)*doseI.DoseGridScaling,[180,180],'LineColor','r','LineWidth',3,'ShowText','on');set(gca, 'YDir','reverse');
        [area,~,~] = Contour2Area(contour(dose(:,:,i)*doseI.DoseGridScaling,[180,180]));
        title(['area: ',num2str(sum(area,'all'))])
%         total_area_inter = total_area_inter + area;
        col = [col;[i, sum(area,'all')]];
    end
    
%     waitforbuttonpress
    close
end 

for i = 1:length(col)-1
    volu_inter = (col(i,2) + col(i+1,2))*5/2;
    total_area_inter = total_area_inter + volu_inter;
end 
disp(sum(total_area_inter/1000,'all'))
%%
fprintf(['P1 prostate volume diff = ',num2str(sum(total_area_inter/1000,'all')-sum(total_area_pre/1000,'all')),'ml\n']);

%% P2 pre
dose = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat2_Preplan\DO001.dcm'));
doseI = dicominfo('C:\Users\SA\Desktop\0531\Thesis\Data\Pat2_Preplan\DO001.dcm');
total_area_pre = 0;
col = [];
for i = 1:size(dose,3)
    figure;
    imshow(dose(:,:,i),[]);axis on; hold on
    contour(dose(:,:,i)*doseI.DoseGridScaling,'ShowText','on');set(gca, 'YDir','reverse');
    if sum((dose(:,:,i)*doseI.DoseGridScaling)>=175,'all') >=1
        hold on
        contour(dose(:,:,i)*doseI.DoseGridScaling,[175,175],'LineColor','r','LineWidth',3,'ShowText','on');set(gca, 'YDir','reverse');
        [area,~,~] = Contour2Area(contour(dose(:,:,i)*doseI.DoseGridScaling,[175,175]));
        title(['area: ',num2str(area)])
%         total_area_pre = total_area_pre + area;
        col = [col;[i, sum(area,'all')]];
    end
%     waitforbuttonpress
    close
end 

for i = 1:length(col)-1
    volu_inter = (col(i,2) + col(i+1,2))*5/2;
    total_area_pre = total_area_pre + volu_inter;
end 
disp(sum(total_area_pre/1000,'all'))

%% P2 interop
dose = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat2_IntraoperativeTrementPlan\DO001.dcm'));
doseI = dicominfo('C:\Users\SA\Desktop\0531\Thesis\Data\Pat2_IntraoperativeTrementPlan\DO001.dcm');
total_area_inter = 0;
col = [];
for i = 1:size(dose,3)
    figure;
    imshow(dose(:,:,i),[]);axis on; hold on
    contour(dose(:,:,i)*doseI.DoseGridScaling,'ShowText','on');set(gca, 'YDir','reverse');
    if sum((dose(:,:,i)*doseI.DoseGridScaling)>=180,'all') >=1
        hold on
        contour(dose(:,:,i)*doseI.DoseGridScaling,[180,180],'LineColor','r','LineWidth',3,'ShowText','on');set(gca, 'YDir','reverse');
        [area,~,~] = Contour2Area(contour(dose(:,:,i)*doseI.DoseGridScaling,[180,180]));
        title(['area: ',num2str(sum(area,'all'))])
%         total_area_inter = total_area_inter + area;
        col = [col;[i, sum(area,'all')]];
    end
    
%     waitforbuttonpress
    close
end 

for i = 1:length(col)-1
    volu_inter = (col(i,2) + col(i+1,2))*5/2;
    total_area_inter = total_area_inter + volu_inter;
end 
disp(sum(total_area_inter/1000,'all'))
%%
fprintf(['P2 prostate volume diff = ',num2str(sum(total_area_inter/1000,'all')-sum(total_area_pre/1000,'all')),'ml\n']);