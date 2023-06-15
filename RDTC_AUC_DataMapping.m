%Import 4D data and crop prior to this step. 
datafc2=double(data);

% -------------------------------------------------------------------------
% Derive functional maps 
% -------------------------------------------------------------------------

%Draw regions of interest over the entire kidney and the medullary/renal
%pelvis (MRP) areas
AUCLobe=[];
SlopeLobe = [];
R2sLobe = [];

for j=1:size(datawfc2,3)
    title('Draw whole left kidney')
    [wlroi(:,:,j)]=roi_2d(datawfc2(:,:,j,5), 25000);
    title('Draw whole middle left kidney')
    [mlroi(:,:,j)]=roi_2d(datawfc2(:,:,j,5), 25000);
    title('Draw whole right kidney')
    [wrroi(:,:,j)]=roi_2d(datawfc2(:,:,j,5), 25000);
    title('Dra whole right middle kidney')
    [mrroi(:,:,j)]=roi_2d(datawfc2(:,:,j,5), 25000);
end

%subtracts the logical ROI matrices to determine the left and right
%cortex/MRP regions 
for j=1:size(datawfc2,3)
LobeROI(:,:,j) = logical((wlroi(:,:,j)-mlroi(:,:,j)) +(wrroi(:,:,j)-mrroi(:,:,j)));
MedullaROI(:,:,j) = logical(mlroi(:,:,j)+mrroi(:,:,j));
LeftLobeROI(:,:,j) = logical(wlroi(:,:,j)-mlroi(:,:,j));
LeftMedullaROI(:,:,j) = logical(mlroi(:,:,j));
RightLobeROI(:,:,j) = logical(wrroi(:,:,j)-mrroi(:,:,j));
RightMedullaROI(:,:,j) = logical(mrroi(:,:,j));
FullROI(:,:,j) = logical(wlroi(:,:,j)+wrroi(:,:,j));
LeftROI(:,:,j) = logical(wlroi(:,:,j));
RightROI(:,:,j) = logical(wrroi(:,:,j));
end

%Sets the time vector for later. 2.5 is minutes of each scan interval so
%this can be changed for whatever scan time you're runnning. 
time(1) = 0;
for j = 2:nscans
    time(j) = 2.5 + time(j-1);
end

%This is for the AUC dataset baselining, it subtracts each pixel at
%each timepoint in each slice from that of t=0 (slice 1) 
for k = 1:size(datawfc2,1)
    for l=1:size(datawfc2,2)
        for m=1:size(datawfc2,3)
            datawfcS(k,l,m,:) = (datawfc2(k,l,m,:)) - (datawfc2(k,l,m,1)); 
        end
    end
end

%zeroes some of the negative values in the AUC if they occur. You can't
%have a negative change in contrast from baseline. 
datawfcS(datawfcS<=0) = 0;

%Gives the AUC value per pixel from the transformed data. 
for k=1:size(datawfc2,1)
    for l=1:size(datawfc2,2)
        for m=1:size(datawfc2,3)
           mapAUC(k,l,m) = trapz(datawfcS(k,l,m,:),4);
        end
    end
end

%Renal decay time constant (not using the AUC dataset, using the original
%starting data set). This performs the linear regression from the highest
%point of contrast (t=2) and the time point where baseline is restored (17
%scans in this case) by taking the natural log of the y data (intensity)
%and the normal x data (time vector from earlier). This outputs R^2 values,
%RDTC values, and CI values per pixel for the ROI on each slice. 
for k=1:size(datawfc2,1)
    disp(['Looping ',num2str(k),' of total ', num2str(size(datawfc2,1))])  
    for l=1:size(datawfc2,2)
        for m=1:size(datawfc2,3)
            xdata=squeeze(time(2:17))';                         
            ydata=squeeze(log(datawfc2(k,l,m,2:17)));     
            [linCoeff, Sw, mu] = polyfit(xdata, ydata, 1);  
            [yfit,delta] = polyval(linCoeff, xdata, Sw, mu);
            SOStot = sum((ydata-mean(ydata)).^2);          
            SOSres = sum((ydata-yfit).^2);               
            Rsquare(k,l,m) = 1-SOSres/SOStot;        
            CI95(k,l,m,:)=2.*delta;                       
            mapslope(k,l,m)=linCoeff(1);                
            mapy(k,l,m)=linCoeff(2);                     
            clear xdata ydata linCoeff Sw mu yfit delta SOStot SOSres CI95 
        end
    end
end

%This makes any values of 0 in R2 and slope (RDTC) a value extremely close
%to 0. This is because when this is mapped later, matlab will read 0 as no
%data and leave a black/clear spot on the image when we try to overlay.
%This prevents that without compromising the numerical data quality. 
Rsquare(isnan(Rsquare))=0.0000001;
mapslope(isnan(mapslope))=0.0000001;

%This is the data writing portion from the transformed data. 
AUCLobe=[];
SlopeLobe = [];
R2sLobe = [];
AUCMedulla=[];
SlopeMedulla = [];
R2sMedulla = [];

AUCLLobe=[];
SlopeLLobe = [];
R2sLLobe = [];
AUCLMedulla=[];
SlopeLMedulla = [];
R2sLMedulla = [];

AUCRLobe=[];
SlopeRLobe = [];
R2sRLobe = [];
AUCRMedulla=[];
SlopeRMedulla = [];
R2sRMedulla = [];

AUCFull=[];
SlopeFull = [];
R2sFull = [];
AUCLeft=[];
SlopeLeft = [];
R2sLeft = [];
AUCRight=[];
SlopeRight = [];
R2sRight = [];

%applies the ROI overtop of the whole transformed data matrix just to give
%the numerical data for the ROI. 
for j=1:size(datawfc2,3)
    R2sLobe = [R2sLobe;nonzeros(squeeze(LobeROI(:,:,j).*Rsquare(:,:,j)))];
    R2sMedulla = [R2sMedulla;nonzeros(squeeze(MedullaROI(:,:,j).*Rsquare(:,:,j)))];
    SlopeLobe = [SlopeLobe;nonzeros(squeeze(LobeROI(:,:,j).*mapslope(:,:,j)))];
    SlopeMedulla = [SlopeMedulla;nonzeros(squeeze(MedullaROI(:,:,j).*mapslope(:,:,j)))];
    AUCLobe = [AUCLobe;nonzeros(squeeze(LobeROI(:,:,j).*mapAUC(:,:,j)))];
    AUCMedulla = [AUCMedulla;nonzeros(squeeze(MedullaROI(:,:,j).*mapAUC(:,:,j)))];
end

for j=1:size(datawfc2,3)
    R2sLLobe = [R2sLLobe;nonzeros(squeeze(LeftLobeROI(:,:,j).*Rsquare(:,:,j)))];
    R2sLMedulla = [R2sLMedulla;nonzeros(squeeze(LeftMedullaROI(:,:,j).*Rsquare(:,:,j)))];
    SlopeLLobe = [SlopeLLobe;nonzeros(squeeze(LeftLobeROI(:,:,j).*mapslope(:,:,j)))];
    SlopeLMedulla = [SlopeLMedulla;nonzeros(squeeze(LeftMedullaROI(:,:,j).*mapslope(:,:,j)))];
    AUCLLobe = [AUCLLobe;nonzeros(squeeze(LeftLobeROI(:,:,j).*mapAUC(:,:,j)))];
    AUCLMedulla = [AUCLMedulla;nonzeros(squeeze(LeftMedullaROI(:,:,j).*mapAUC(:,:,j)))];
end

for j=1:size(datawfc2,3)
    R2sRLobe = [R2sRLobe;nonzeros(squeeze(RightLobeROI(:,:,j).*Rsquare(:,:,j)))];
    R2sRMedulla = [R2sRMedulla;nonzeros(squeeze(RightMedullaROI(:,:,j).*Rsquare(:,:,j)))];
    SlopeRLobe = [SlopeRLobe;nonzeros(squeeze(RightLobeROI(:,:,j).*mapslope(:,:,j)))];
    SlopeRMedulla = [SlopeRMedulla;nonzeros(squeeze(RightMedullaROI(:,:,j).*mapslope(:,:,j)))];
    AUCRLobe = [AUCRLobe;nonzeros(squeeze(RightLobeROI(:,:,j).*mapAUC(:,:,j)))];
    AUCRMedulla = [AUCRMedulla;nonzeros(squeeze(RightMedullaROI(:,:,j).*mapAUC(:,:,j)))];
end

for j=1:size(datawfc2,3)
    R2sFull = [R2sFull;nonzeros(squeeze(FullROI(:,:,j).*Rsquare(:,:,j)))];
    R2sLeft = [R2sLeft;nonzeros(squeeze(LeftROI(:,:,j).*Rsquare(:,:,j)))];
    R2sRight = [R2sRight;nonzeros(squeeze(RightROI(:,:,j).*Rsquare(:,:,j)))];
    SlopeFull = [SlopeFull;nonzeros(squeeze(FullROI(:,:,j).*mapslope(:,:,j)))];
    SlopeLeft = [SlopeLeft;nonzeros(squeeze(LeftROI(:,:,j).*mapslope(:,:,j)))];
    SlopeRight = [SlopeRight;nonzeros(squeeze(RightROI(:,:,j).*mapslope(:,:,j)))];
    AUCFull = [AUCFull;nonzeros(squeeze(FullROI(:,:,j).*mapAUC(:,:,j)))];
    AUCLeft = [AUCLeft;nonzeros(squeeze(LeftROI(:,:,j).*mapAUC(:,:,j)))];
    AUCRight = [AUCRight;nonzeros(squeeze(RightROI(:,:,j).*mapAUC(:,:,j)))];
end

%Writes all of these data matrices to excel. 
writematrix(AUCLobe, 'AUCLobe.csv', 'WriteMode', 'append');
writematrix(SlopeLobe, 'SlopeLobe.csv', 'WriteMode', 'append');
writematrix(R2sLobe, 'R2sLobe.csv', 'WriteMode', 'append');
writematrix(AUCMedulla, 'AUCMedulla.csv', 'WriteMode', 'append');
writematrix(SlopeMedulla, 'SlopeMedulla.csv', 'WriteMode', 'append');
writematrix(R2sMedulla, 'R2sMedulla.csv', 'WriteMode', 'append');

writematrix(AUCLLobe, 'L-AUCLobe.csv', 'WriteMode', 'append');
writematrix(SlopeLLobe, 'L-SlopeLobe.csv', 'WriteMode', 'append');
writematrix(R2sLLobe, 'L-R2sLobe.csv', 'WriteMode', 'append');
writematrix(AUCLMedulla, 'L-AUCMedulla.csv', 'WriteMode', 'append');
writematrix(SlopeLMedulla, 'L-SlopeMedulla.csv', 'WriteMode', 'append');
writematrix(R2sLMedulla, 'L-R2sMedulla.csv', 'WriteMode', 'append');

writematrix(AUCRLobe, 'R-AUCLobe.csv', 'WriteMode', 'append');
writematrix(SlopeRLobe, 'R-SlopeLobe.csv', 'WriteMode', 'append');
writematrix(R2sRLobe, 'R-R2sLobe.csv', 'WriteMode', 'append');
writematrix(AUCRMedulla, 'R-AUCMedulla.csv', 'WriteMode', 'append');
writematrix(SlopeRMedulla, 'R-SlopeMedulla.csv', 'WriteMode', 'append');
writematrix(R2sRMedulla, 'R-R2sMedulla.csv', 'WriteMode', 'append');

writematrix(AUCFull, 'AUCFull.csv', 'WriteMode', 'append');
writematrix(SlopeFull, 'SlopeFull.csv', 'WriteMode', 'append');
writematrix(R2sFull, 'R2sFull.csv', 'WriteMode', 'append');
writematrix(AUCLeft, 'AUCLeft.csv', 'WriteMode', 'append');
writematrix(SlopeLeft, 'SlopeLeft.csv', 'WriteMode', 'append');
writematrix(R2sLeft, 'R2sLeft.csv', 'WriteMode', 'append');
writematrix(AUCRight, 'AUCRight.csv', 'WriteMode', 'append');
writematrix(SlopeRight, 'SlopeRight.csv', 'WriteMode', 'append');
writematrix(R2sRight, 'R2sRight.csv', 'WriteMode', 'append');

%This gives the average intensity values for the untransformed images at
%each time point. Writes to a matrix. 
AvgIntensityLobe=[];
AvgIntensityMed=[];
AvgIntensity=[];

for k = 1:size(datawfc2,3)
    for j=1:size(datawfc2,4)
    AvgIntensityLobe(j,k+1) = squeeze(mean(nonzeros(LobeROI(:,:,k).*datawfc2(:,:,k,j))));
    AvgIntensityMed(j,k+1) = squeeze(mean(nonzeros(MedullaROI(:,:,k).*datawfc2(:,:,k,j))));
    AvgIntensityFull(j,k+1) = squeeze(mean(nonzeros(FullROI(:,:,k).*datawfc2(:,:,k,j))));
    end
end

for j = 1:size(datawfc2,4)
    AvgIntensityMed(j,size(datawfc2,3)+2) = sum(AvgIntensityMed(j,2:1+(size(datawfc2,3))));
    AvgIntensityLobe(j,size(datawfc2,3)+2) = sum(AvgIntensityLobe(j,2:1+(size(datawfc2,3))));
    AvgIntensityFull(j,size(datawfc2,3)+2) = sum(AvgIntensityFull(j,2:1+(size(datawfc2,3))));
end 

AvgIntensityMed(1) = 0;
for j = 2:nscans
    AvgIntensityMed(j,1) = 2.5 + AvgIntensityMed(j-1);
end

AvgIntensityLobe(1) = 0;
for j = 2:nscans
    AvgIntensityLobe(j,1) = 2.5 + AvgIntensityLobe(j-1);
end

AvgIntensityFull(1) = 0;
for j = 2:nscans
    AvgIntensityFull(j,1) = 2.5 + AvgIntensityFull(j-1);
end

writematrix(AvgIntensityMed, 'AvgIntensityMed.csv', 'WriteMode', 'append');
writematrix(AvgIntensityLobe,  'AvgIntensityLobe.csv', 'WriteMode', 'append');
writematrix(AvgIntensityFull, 'AvgIntensityFull.csv', 'WriteMode', 'append');

%Making the overlayed map images from the written transformed matrices.
%Overlays the mapped (AUC or RDTC) ROI on top of the original grayscale MRI
%image. Also saves the images as matlab files for later saving as tiff. 
cropslope = mapslope;
cropauc = mapAUC;
cropslope(cropslope==0) = NaN;
cropauc(cropauc<=0) = NaN;

for j=1:size(datawfc2,3)
    h=figure(j);
    imagesc(mapslope(:,:,j), [-0.3 0.3]), colormap turbo;
    colorbar;
    saveas(h,sprintf('Slope%d.fig',j));
    reset(j);
    close(j);
end

for j=1:size(datawfc2,3)
    h=figure(j);
    imagesc(mapAUC(:,:,j), [0 250000]), colormap jet
    colorbar
    saveas(h,sprintf('AUC%d.fig',j));
    reset(j);
    close(j);
end

for j=1:size(data,3)
    y = size(data(:,:,j),1);
    x = size(data(:,:,j),2);
end

for j=1:size(datawfc2,3)
    hf = figure(j);
    h1 = axes;colormap(h1,'gray');
    p1=imagesc(x,y,datawfc2(:,:,j))
    set(h1,'xdir','normal');
    colormap(h1, 'gray');
    % Foreground image
    h2=axes;
    set(h2,'ydir','normal');
    p2=imagesc(x,y,cropslope(:,:,j),"AlphaData", 0.40);
    colormap(h2, 'turbo')
    set(h2,'color','none','visible','off')
    clim([-0.3 0.3]);
    linkaxes([h1 h2])
    saveas(hf,sprintf('SlopeFig%d.tiff',j));
    reset(hf)
    close(hf);
end

for j=1:size(datawfc2,3)
    hf = figure(j+50);
    h1 = axes;colormap(h1,'gray');
    p1=imagesc(x,y,datawfc2(:,:,j))
    set(h1,'xdir','normal');
    colormap(h1, 'gray');
    % Foreground image
    h2=axes;
    set(h2,'ydir','normal');
    p2=imagesc(x,y,cropauc(:,:,j),"AlphaData", 0.40);
    colormap(h2, 'jet')
    set(h2,'color','none','visible','off')
    clim([0 250000]);
    linkaxes([h1 h2])
    saveas(hf,sprintf('AUCFig%d.tiff',j));
    reset(hf)
    close(hf);
end

for j=1:size(datawfc2,3)
    hf = figure(j);
    h1 = axes;colormap(h1,'gray');
    p1=imagesc(x,y,datawfc2(:,:,j,2))
    set(h1,'xdir','normal');
    colormap(h1, 'gray');
    clim([0 30000]);
    saveas(hf,sprintf('GSFig%d.tiff',j));
    reset(hf)
    close(hf);
end

