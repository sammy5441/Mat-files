% % <summary>
% % Task 3.0: TG-43 Comparision
% % </summary>
% % <remarks>
% %  Author:            SG
% %                     (C) Heidelberg University
% %  Project name:      Master Thesis : Seed train positions (orientation and trajectories) of
% %                     the pre-plan and the final intra-operative implant using the dicom files:
% %
% %  Date:              2021-05-11
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

close all;
% clear all;
clc;

%% Source Specifications for BrachySource? 125I Seeds
%% Radial Dose Function
%% Point Source Approximation


x = [0.1 0.15 0.25 0.5 0.75 1.00 1.5 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0];
y = [0.544 0.7 0.876 0.999 1.013 1.00 0.943 0.864 0.698 0.546 0.420 0.318 0.239 0.178 0.133 0.0980];


figure
plot(x,y,'--')
grid on

%Seed Of interest 1
dose = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_IntraoperativeTrementPlan\DO001.dcm'));
doseI = dicominfo('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_IntraoperativeTrementPlan\DO001.dcm');
doseRes = 0.122751/doseI.DoseGridScaling;
yD = flip(double(dose(43,15:25,4))); % Coordinates depend on the Dose Map, must be chosen manualy
xD = doseRes/2:doseRes:doseRes*(length(yD));

hold on
plot(xD/100,yD/max(yD))

%Seed Of interest 2

yD2 = (double(dose(46,67:77,4))); % Coordinates depend on the Dose Map, must be chosen manualy
xD2 = doseRes/2:doseRes:doseRes*(length(yD2));

plot(xD2/100,yD2/max(yD2))


%Seed Of interest 3 (MEAN)
yDm = mean([yD/max(yD);yD2/max(yD2)],1); % Coordinates depend on the Dose Map, must be chosen manualy

plot(xD/100,yDm/max(yDm))


xlabel('Depth (cm)')
ylabel('g(r)')
title('Radial dose functions')
legend('Literature values','Seed 1','Seed 2','Mean')


