clc,close all, clear all
%%
D = "C:\Users\SA\Desktop\0531\Thesis\Data\Mat files\dicom-dict-iotp.txt"
seedinfo = dicominfo('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_IntraoperativeTrementPlan/PL001.dcm','dictionary',D);
seeds = struct2array(seedinfo.ApplicationSetupSequence);

%%
dose = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_IntraoperativeTrementPlan/DO001.dcm'));
for i=1:14
    figure
    img = dose(:,:,i);
    subplot(1,2,1)
    imshow(img,[]);
    img_16 = uint16((img-min(img,[],'all'))/(max(img,[],'all')-min(img,[],'all'))*65535);
    [mask,center] = imsegkmeans(img_16,3);
    subplot(1,2,2)
    imshow(mask,[]);
    waitforbuttonpress
    close
end 
%% how to visulize the SEED from dose map
% just stack 14 slice to 1
dosemap = sum(dose,3);
subplot(1,3,1)
imshow(dosemap,[]);
ma = max(dosemap,[],'all');
mi = min(dosemap,[],'all');
dosemap_16 = uint16((dosemap-mi)/(ma-mi)*65536);
[L,center_] = imsegkmeans(dosemap_16,3);
subplot(1,3,2)
imshow(L,[])
subplot(1,3,3)
imshow(dosemap.*(L==2),[])
%% intensity map
% could 

spot = {};
for i = 1:80
    [~, loc] = findpeaks([single(dosemap_16(i,:))]);
    disp(size(loc))
    if size(loc)==[1,1]
%         spot{i} = [];
        spot{i} = loc;
    else
        spot{i} = loc;
    end
end

spot_v = {};
for i = 1:80
    [~, loc] = findpeaks([single(dosemap_16(:,i))]);
    disp(size(loc))
    if size(loc)==[1,1]
%         spot{i} = [];
        spot_v{i} = loc;
    else
        spot_v{i} = loc;
    end
end


testmap = zeros(80,80);
for i = 1:80
    testmap(i,spot{i}) = testmap(i,spot{i}) + 1;
end 
for i = 1:80
    testmap(spot_v{i},i) = testmap(spot_v{i},i) + 1;
end 
figure
subplot(1,3,1)
imshow(dosemap_16,[])
subplot(1,3,2)
imshow(testmap,[])
subplot(1,3,3)
imshow(testmap==2,[])

% [r,c] = find(testmap==2);
% y_all = []
% figure
% for i = 1:length(r)
%     yD = flip(double(dose(r_t(i),c_t(i)-4:c_t(i)+4,z_t(i)))); % Coordinates depend on the Dose Map, must be chosen manualy
%     xD = doseRes/2:doseRes:doseRes*(length(yD));
%     hold on
%     plot(xD/100,yD/max(yD))
%     y_all = [y_all;yD/max(yD)];
% end 
% yDm = mean(y_all,1);
% plot(xD/100,yDm/max(yDm))
% hold on
% x = [0.1 0.15 0.25 0.5 0.75 1.00 1.5 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0];
% y = [0.544 0.7 0.876 0.999 1.013 1.00 0.943 0.864 0.698 0.546 0.420 0.318 0.239 0.178 0.133 0.0980];
% plot(x,y,'--')
% legend
%%
r_t=[]; c_t=[]; z_t=[];
totalmap = zeros(80,80,14);
for z_slice = 3:12
    dosemap_16 = rot90(rot90(dose(:,:,z_slice)));
    spot = {};
    for i = 1:80
        [~, loc] = findpeaks([single(dosemap_16(i,:))]);
%         disp(size(loc))
        if size(loc)==[1,1]
    %         spot{i} = [];
            spot{i} = loc;
        else
            spot{i} = loc;
        end
    end

    spot_v = {};
    for i = 1:80
        [~, loc] = findpeaks([single(dosemap_16(:,i))]);
        disp(size(loc))
        if size(loc)==[1,1]
    %         spot{i} = [];
            spot_v{i} = loc;
        else
            spot_v{i} = loc;
        end
    end


    testmap = zeros(80,80);
    for i = 1:80
        testmap(i,spot{i}) = testmap(i,spot{i}) + 1;
    end 
    for i = 1:80
        testmap(spot_v{i},i) = testmap(spot_v{i},i) + 1;
    end 
%     figure
%     subplot(1,3,1)
%     imshow(dosemap_16,[])
%     subplot(1,3,2)
%     imshow(testmap,[])
%     subplot(1,3,3)
%     imshow(testmap==2,[])
%     totalmap(:,:,z_slice) = testmap==2;
    [r,c] = find(testmap==2);
    z = ones(1,length(r))*z_slice;
    r_t = [r_t,rot90(r)];
    c_t = [c_t,rot90(c)];
    z_t = [z_t,z];
    
end 
figure
scatter3(c_t, r_t, z_t, 'g');

%%

dose = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_IntraoperativeTrementPlan\DO001.dcm'));
doseI = dicominfo('C:\Users\SA\Desktop\0531\Thesis\Data\Pat1_IntraoperativeTrementPlan\DO001.dcm','dictionary',D);
doseRes = 0.122751/doseI.DoseGridScaling;
y_all = [];
for i = 1:131
    yD = flip(double(dose(r_t(i),c_t(i)-5:c_t(i)+5,z_t(i)))); % Coordinates depend on the Dose Map, must be chosen manualy
    xD = doseRes/2:doseRes:doseRes*(length(yD));
    hold on
    plot(xD/100,yD/max(yD))
    y_all = [y_all;yD/max(yD)];
end 
yDm = mean(y_all,1);
figure
plot(xD/100,yDm/max(yDm))
hold on 
x = [0.1 0.15 0.25 0.5 0.75 1.00 1.5 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0];
y = [0.544 0.7 0.876 0.999 1.013 1.00 0.943 0.864 0.698 0.546 0.420 0.318 0.239 0.178 0.133 0.0980];
plot(x,y,'--')
grid on
%%
img = dicomread('./MR008.dcm');
imshow(img,[]);