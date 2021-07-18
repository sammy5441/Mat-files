clc, clear all, close all
%%
dosePre = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat2_Preplan\DO001.dcm'));
dosePost = squeeze(dicomread('C:\Users\SA\Desktop\0531\Thesis\Data\Pat2_IntraoperativeTrementPlan\DO001.dcm'));
doseDelta = abs(dosePost-dosePre);


for i = 1:10
    figure;
    subplot(1,3,1)
    imshow(dosePre(:,:,i),[]);axis on; hold on
    colormap(jet(256));
    contour(dosePre(:,:,i)*doseI.DoseGridScaling,'ShowText','on');set(gca, 'YDir','reverse');
    subplot(1,3,2)
    imshow(dosePost(:,:,i),[]);axis on; hold on
    colormap(jet(256));
    contour(dosePost(:,:,i)*doseI.DoseGridScaling,'ShowText','on');set(gca, 'YDir','reverse');
    subplot(1,3,3)
    imshow(doseDelta(:,:,i),[]);axis on; hold on
    colormap(jet(256));
    contour(doseDelta(:,:,i)*doseI.DoseGridScaling,'ShowText','on');set(gca, 'YDir','reverse');
    i=i+3;%k=k+3;
end

doseDeltaInterp = interp3(cast(doseDelta,"double"));
figure
a = slice(doseDeltaInterp,[],[],[1:1:10]);
shading flat
alpha(a, 0.5)

figure
dose = {dosePre,dosePost, doseDelta};
for  l = 1:1:3
    doseDeltaInterp = interp3(cast(dose{l},"double"));
    subplot(1,3,l)
    aval = 0.0;
    i = 1;
    col = ["white",'blue',"green","yellow","red"];
    for j = 0:max(doseDeltaInterp,[],"all")/5:max(doseDeltaInterp,[],"all")*0.99
        b = patch(isosurface(doseDeltaInterp,j),"FaceAlpha",aval);
        
        b.FaceColor = col(i);
        b.EdgeColor = 'none';
        daspect([1 1 1])
        view(3);
        axis tight
        camlight
        lighting gouraud
        grid minor
        i = i +1;
        aval = aval + 0.2;
    end
end

%Pull Quad Data
quadStats = cell(4,4);
q = 1;
figure
for i = [1,40]
    for j = [1,40]
        subplot(2,2,q)
        quad = cast(doseDelta(i:i+40,j:j+40,:),"double");
        hist(reshape(quad,[1,41*41*10]))
        quadStats{q,1} = quad;
        quadStats{q,3} =std(quad,[],"all");
        quadStats{q,2} =mean(quad,"all");
        quadStats{q,4} =median(quad,"all");
        q = q+1;
    end
end