% % <summary>
% % Task 3.0: DVH Comparison
% % </summary>
% % <remarks>
% %  Author:            SG
% %                     (C) Heidelberg University
% %  Project name:      Master Thesis : Seed train positions (orientation and trajectories) of
% %                     the pre-plan and the final intra-operative implant using the dicom files:
% %
% %  Date:              2021-04-21
% % </remarks>
% % % % % % % % % % % % % % % % % % % %


%% file
% switch lower( getenv( 'COMPUTERNAME' ) )
%     case 'computername' % your computer name
%         path_proj = 'C:\----'; % your folder location
%     case 'desktop-ae86p3r'
%         path_proj = 'C:\Users\Samson\Documents\Data';
%     otherwise
%         return
% end



clc, close all

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

% yr_plan = [yr(1:end)];
% Vol_yr_plan = trapz(yr)
% cvol_yr_plan = cumtrapz(yr);
% 
% yr_io = [yr_io(1:end)];
% Vol_yr_io = trapz(yr_io)
% cvol_yr_io = cumtrapz(yr_io);
% 
% yu_plan = [yu(1:end)];
% Vol_yu_plan = trapz(yu)
% cvol_yu_plan = cumtrapz(yu);
% 
% yu_io = [yu_io(1:end)];
% Vol_yu_io = trapz(yu_io)
% cvol_yu_io = cumtrapz(yu_io);

subplot(1,2,2)
plot(cvol_plan_p1,'r')
hold on
plot(cvol_io_p1,'r--')
hold on
plot(cvol_plan_p2,'b')
hold on
plot(cvol_io_p2,'b--')
% hold on
% plot(cvol_yr_plan)
% hold on
% plot(cvol_yr_io)
% hold on
% plot(cvol_yu_plan)
% hold on
% plot(cvol_yu_io)
legend("Prostate Plan Vol_ P1","Prostate IntraOp Vol_ P1","Prostate Plan Vol_ P2","Prostate IntraOp Vol_ P2",'Location','southeast')
title('Cumulative vol')
xlabel('Dose (Gy)')
ylabel('Vol')

% f = fit(x1,y1,'poly9')
% plot(f,x,y)
