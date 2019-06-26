% This code is used to obtain infrastructure information for holes and 
% particles in tomography data.
% in Avizo, you should convert image type to 32-bit signed and export 
% data as 3D raw images. the image stack size should be recorded below.
%
% @Author : hwchai_PIMS    04:05 PM   06/24/2019

%[M,N,K] = deal(300,317,247);  % image stack size

PixelSize = 0.87;  % pixel size 0.87 um

% filter_coor = {'volume > 1000','barycenterZ > 100'};

%% -------------------------------------------------------------------------
[I1,filename] = import_raw();  % imread 3D image stack

[vol,gray_value] = separate_grains(I1);
% vol is coordination of each pixel_j in each partical_i  vol{i,1}(1:3,j)
% gray_value(i,2) save the gray value of each partical_i in original image

structu = gyration_tensor(vol,gray_value,PixelSize);
% The structure parameters of each partical will be saved in 'structu'

histo = histogram_(structu);
% Statistics for structural parameters

visual_(histo,filename);
% visualization statistics

clear gray_value M N K PixelSize vol I1 filename;

%% -------------------------------------------------------------------------
function [ histo ] = histogram_(structu)

histo = struct('Volume3d',[],'EqDiameter',[],'Sphericity',[],...
	'ElongationIndex',[],'FlatnessIndex',[],'Convexity',[]);

for i = 1:size(structu,2),
Volume3d(i,1) = structu(i).Volume3d;
EqDiameter(i,1) = structu(i).EqDiameter;
Sphericity(i,1) = structu(i).Sphericity;
ElongationIndex(i,1) = structu(i).ElongationIndex;
FlatnessIndex(i,1) = structu(i).FlatnessIndex;
Convexity(i,1) = structu(i).Convexity;
end

edge_Volume3d = 0:max(Volume3d)/20:max(Volume3d);
% statistical interval of Volume3d
histo.Volume3d(:,1) = (edge_Volume3d(1:end-1)+edge_Volume3d(2:end))/2;
h = histogram(Volume3d,edge_Volume3d);
histo.Volume3d(:,2) = h.Values/sum(h.Values);
histo.Volume3d(:,3) = cumsum(histo.Volume3d(:,2));


edge_EqDiameter = 0:max(EqDiameter)/20:max(EqDiameter);
% statistical interval of EqDiameter
histo.EqDiameter(:,1) = (edge_EqDiameter(1:end-1)+edge_EqDiameter(2:end))/2;
h = histogram(EqDiameter,edge_EqDiameter);
histo.EqDiameter(:,2) = h.Values/sum(h.Values);
histo.EqDiameter(:,3) = cumsum(histo.EqDiameter(:,2));


edge_structu = 0:0.05:1;
% statistical interval of S,EI,FI,Cx
histo.Sphericity(:,1) = (edge_structu(1:end-1)+edge_structu(2:end))/2;
histo.ElongationIndex(:,1) = (edge_structu(1:end-1)+edge_structu(2:end))/2;
histo.FlatnessIndex(:,1) = (edge_structu(1:end-1)+edge_structu(2:end))/2;
histo.Convexity(:,1) = (edge_structu(1:end-1)+edge_structu(2:end))/2;

h = histogram(Sphericity,edge_structu);
histo.Sphericity(:,2) = h.Values/sum(h.Values);
histo.Sphericity(:,3) = cumsum(histo.Sphericity(:,2));

h = histogram(ElongationIndex,edge_structu);
histo.ElongationIndex(:,2) = h.Values/sum(h.Values);
histo.ElongationIndex(:,3) = cumsum(histo.ElongationIndex(:,2));

h = histogram(FlatnessIndex,edge_structu);
histo.FlatnessIndex(:,2) = h.Values/sum(h.Values);
histo.FlatnessIndex(:,3) = cumsum(histo.FlatnessIndex(:,2));

h = histogram(Convexity,edge_structu);
histo.Convexity(:,2) = h.Values/sum(h.Values);
histo.Convexity(:,3) = cumsum(histo.Convexity(:,2));

close all

end

%% -------------------------------------------------------------------------
function [ ] = visual_(histo,filename)

fig = figure; set(fig,'units','normalized','position',[0.2 0.2 0.6 0.6]);

subplot(2,3,1);
[ax,p1,p2] = plotyy(histo.Volume3d(:,1),histo.Volume3d(:,2),...
	histo.Volume3d(:,1),histo.Volume3d(:,3));
set(ax,'Linewidth',1);
xlabel(ax(1),'Volume (um^3)','fontsize',10);
ylabel(ax(1),'Frequency','fontsize',10);
ylabel(ax(2),'Cumulative probability','fontsize',10);
set(ax(1),'FontName','Helvetica','FontSize',8);
set(ax(2),'FontName','Helvetica','FontSize',8);
grid on;

subplot(2,3,2);
[ax,p1,p2] = plotyy(histo.Sphericity(:,1),histo.Sphericity(:,2),...
	histo.Sphericity(:,1),histo.Sphericity(:,3));
set(ax,'Linewidth',1);
xlabel(ax(1),'Sphericity','fontsize',10);
ylabel(ax(1),'Frequency','fontsize',10);
ylabel(ax(2),'Cumulative probability','fontsize',10);
set(ax(1),'FontName','Helvetica','FontSize',8);
set(ax(2),'FontName','Helvetica','FontSize',8);
grid on;

subplot(2,3,3);
[ax,p1,p2] = plotyy(histo.ElongationIndex(:,1),histo.ElongationIndex(:,2),...
	histo.ElongationIndex(:,1),histo.ElongationIndex(:,3));
set(ax,'Linewidth',1);
xlabel(ax(1),'Elongation Index','fontsize',10);
ylabel(ax(1),'Frequency','fontsize',10);
ylabel(ax(2),'Cumulative probability','fontsize',10);
set(ax(1),'FontName','Helvetica','FontSize',8);
set(ax(2),'FontName','Helvetica','FontSize',8);
grid on;

subplot(2,3,4);
[ax,p1,p2] = plotyy(histo.EqDiameter(:,1),histo.EqDiameter(:,2),...
	histo.EqDiameter(:,1),histo.EqDiameter(:,3));
set(ax,'Linewidth',1);
xlabel(ax(1),'Equivalent Diameter (um)','fontsize',10);
ylabel(ax(1),'Frequency','fontsize',10);
ylabel(ax(2),'Cumulative probability','fontsize',10);
set(ax(1),'FontName','Helvetica','FontSize',8);
set(ax(2),'FontName','Helvetica','FontSize',8);
grid on;

subplot(2,3,5);
[ax,p1,p2] = plotyy(histo.FlatnessIndex(:,1),histo.FlatnessIndex(:,2),...
	histo.FlatnessIndex(:,1),histo.FlatnessIndex(:,3));
set(ax,'Linewidth',1);
xlabel(ax(1),'Flatness Index','fontsize',10);
ylabel(ax(1),'Frequency','fontsize',10);
ylabel(ax(2),'Cumulative probability','fontsize',10);
set(ax(1),'FontName','Helvetica','FontSize',8);
set(ax(2),'FontName','Helvetica','FontSize',8);
grid on;

subplot(2,3,6);
[ax,p1,p2] = plotyy(histo.Convexity(:,1),histo.Convexity(:,2),...
	histo.Convexity(:,1),histo.Convexity(:,3));
set(ax,'Linewidth',1);
xlabel(ax(1),'Convexity','fontsize',10);
ylabel(ax(1),'Frequency','fontsize',10);
ylabel(ax(2),'Cumulative probability','fontsize',10);
set(ax(1),'FontName','Helvetica','FontSize',8);
set(ax(2),'FontName','Helvetica','FontSize',8);
grid on;

print(fig,['histogram_',filename,'.png'],'-dpng','-r200');

end

%% -------------------------------------------------------------------------
function [I1,filename2] = import_raw()
% import 3D raw image stack which with image stack size [M,N,k]
[filename, pathname] = uigetfile('*.raw', 'import raw image');

NameSplit = strsplit(filename,'_');
[M,N,K] = deal(str2num(NameSplit{1,end-2}),str2num(NameSplit{1,end-1}),...
	str2num(NameSplit{1,end}));
filename2 = NameSplit{1,1};

f1 = fopen([pathname, filename], 'r');
I1 = fread(f1, 'int32'); fclose(f1);
len = length(I1);
if K ~= len/(M*N), error('image stack size error !'); end
I1 = reshape(I1,M,N,K);
end

%% -------------------------------------------------------------------------
function [vol,gray_value] = separate_grains(I1)
% separate each grains to sparse matrix
[M,N,K] = size(I1);
% cell_number = max(max(max(I1)));  % particle number
% vol = cell(cell_number,1);        
% coordinates of pixel_i in particle n vol{n,1}(:,i);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
I1 = I1(:); I1(:,2) = (1:size(I1,1))';
I1(find(I1(:,1) == 0),:) = [];
[I1(:,3),I1(:,4),I1(:,5)] = ind2sub([M,N,K],I1(:,2));
I1 = sortrows(I1,1);
I1(2:end,6) = I1(1:end-1,1); I1(:,6) = I1(:,6) - I1(:,1);
edge_ = find(I1(:,6)~=0); 
edge_(end+1,1) = size(I1,1)+1;  % Add the end index of the last partical
cell_number = size(edge_,1)-1;  % particle number
vol = cell(cell_number,1);
gray_value = (1:cell_number)';
for i = 1:size(edge_,1)-1,
	vol{i,1} = I1(edge_(i,1):edge_(i+1,1)-1,3:5);
	gray_value(i,2) = I1(edge_(i,1),1); 
end

end

function [structu] = gyration_tensor(vol,gray_value,PixelSize)
%% VARIABLES OPTIONS
% OUTPUT NOTES
% -------------------------------------------------------------------------
structu = struct('ParticalIndex',[],...         % partical index
                 'GrayValue',[],...             % gray value in original image
                 'BaryCenterX',[],...           % coordination X of barycenter
                 'BaryCenterY',[],...           % coordination Y of barycenter
                 'BaryCenterZ',[],...           % coordination Z of barycenter
	             'Volume3d',[],...              % digital volume
                 'EqDiameter',[],...            % equivalent diameter
                 'Sphericity',[],...            % Sphericity  @Y.Yao JMS  1-A
                 'ElongationIndex',[],...       % EI  @H.Y.Li Powder Techn
                 'FlatnessIndex',[],...         % FI  @H.Y.Li Powder Techn
                 'Convexity',[],...             % CX  @H.Y.Li Powder Techn
	             'GyrationTensor',[],...         % gyration tensor @Y.Yao JMS G
                 'EigenValA',[],...             % R1 (or lambda1) = a^2/5
                 'EigenValB',[],...             % R2 (or lambda2) = b^2/5
                 'EigenValC',[],...             % R3 (or lambda3) = c^2/5
                 'EigenVecA',[],...             % coordination X of R1 Axis
                 'EigenVecB',[],...             % coordination X of R2 Axis
                 'EigenVecC',[],...             % coordination X of R3 Axis
                 'AxisATheta',[],...            % Theta of R1 Axis
                 'AxisAPhi',[],...              % Phi of R1 Axis
                 'AxisBTheta',[],...            % Theta of R2 Axis
                 'AxisBPhi',[],...              % Phi of R2 Axis
                 'AxisCTheta',[],...            % Theta of R3 Axis
                 'AxisCPhi',[]);                % Phi of R3 Axis
%                 'Area3d',[],...               % Area of partical
%                 'Distance',[],
% Gives the shortest edge to edge distance from the current object to its nearest neighbor.
%%

cell_number = size(vol,1);
for i = 1:cell_number,
structu(i).ParticalIndex = i;
structu(i).GrayValue = gray_value(i,2);
structu(i).BaryCenterX = mean(vol{i,1}(:,1)); % coordination for first pixel is 1, not 0
structu(i).BaryCenterY = mean(vol{i,1}(:,2));
structu(i).BaryCenterZ = mean(vol{i,1}(:,3));
structu(i).Volume3d = size(vol{i,1},1);
structu(i).EqDiameter = 2*(3*structu(i).Volume3d/(4*pi))^(1/3);

for m = 1:3, for n = 1:3,
structu(i).GyrationTensor(m,n) = mean((vol{i,1}(:,m)- mean(vol{i,1}(:,m))).*...
	(vol{i,1}(:,n)- mean(vol{i,1}(:,n)))); end,end

[EigenVec,EigenVal] = eig(structu(i).GyrationTensor);

for m = 1:3, EigenVec(4,m) = EigenVal(m,m); end

EigenVec = sortrows(EigenVec',4,'descend')';

EigenVec(1:3,1) = EigenVec(1:3,1)/(abs(EigenVec(3,1))/EigenVec(3,1));
EigenVec(1:3,2) = EigenVec(1:3,2)/(abs(EigenVec(3,2))/EigenVec(3,2));
EigenVec(1:3,3) = EigenVec(1:3,3)/(abs(EigenVec(3,3))/EigenVec(3,3));

structu(i).EigenValA = EigenVec(4,1);
structu(i).EigenValB = EigenVec(4,2);
structu(i).EigenValC = EigenVec(4,3);

structu(i).EigenVecA = EigenVec(1:3,1);
structu(i).EigenVecB = EigenVec(1:3,2);
structu(i).EigenVecC = EigenVec(1:3,3);

structu(i).AxisATheta = 180*asin(sqrt(sum(EigenVec(1:2,1).^2))/sqrt(sum(EigenVec(1:3,1).^2)))/pi;
structu(i).AxisBTheta = 180*asin(sqrt(sum(EigenVec(1:2,2).^2))/sqrt(sum(EigenVec(1:3,2).^2)))/pi;
structu(i).AxisCTheta = 180*asin(sqrt(sum(EigenVec(1:2,3).^2))/sqrt(sum(EigenVec(1:3,3).^2)))/pi;

structu(i).AxisAPhi = 180*atan(EigenVec(1,1)/EigenVec(2,1))/pi + 180*(1-abs(EigenVec(1,1))/EigenVec(1,1))/2;
structu(i).AxisBPhi = 180*atan(EigenVec(1,2)/EigenVec(2,2))/pi + 180*(1-abs(EigenVec(1,2))/EigenVec(1,2))/2;
structu(i).AxisCPhi = 180*atan(EigenVec(1,3)/EigenVec(2,3))/pi + 180*(1-abs(EigenVec(1,3))/EigenVec(1,3))/2;

structu(i).Sphericity = 3*sum(EigenVec(4,:).*EigenVec(4,[3 1 2]))/(sum(EigenVec(4,:)))^2;
% 3*(R1R2+R1R3+R2R3)/(R1+R2+R3)^2. spherical objects will have the values close to 1.

structu(i).ElongationIndex = structu(i).EigenValB/structu(i).EigenValA;
% The ratio of the medium to the largest eigenvalue of gyration tensor. 
% Elongated objects will have small values close to 0.

structu(i).FlatnessIndex = structu(i).EigenValC/structu(i).EigenValB;
% The ratio of the smallest to the medium eigenvalue of the covariance 
% matrix. Flat objects have small values close to 0.

%[triangulation, ConvexVolume] = convhulln(vol{i,1}(:,1),vol{i,1}(:,2),vol{i,1}(:,3));
[triangulation, ConvexVolume] = convhulln(vol{i,1});
structu(i).Convexity = structu(i).Volume3d/ConvexVolume;
% The ratio of the volume of convex hull to original volume.
% BUG : Convexity of those partical which cropped by field of view will
% larger then 1. It's not reasonable.

% correlation from pixel to um
structu(i).BaryCenterX = structu(i).BaryCenterX*PixelSize; 
structu(i).BaryCenterY = structu(i).BaryCenterY*PixelSize;
structu(i).BaryCenterZ = structu(i).BaryCenterZ*PixelSize;

structu(i).Volume3d = structu(i).Volume3d*PixelSize^3;
structu(i).EqDiameter = structu(i).EqDiameter*PixelSize;

structu(i).GyrationTensor = structu(i).GyrationTensor*PixelSize^2;

structu(i).EigenValA = structu(i).EigenValA*PixelSize^2;
structu(i).EigenValB = structu(i).EigenValB*PixelSize^2;
structu(i).EigenValC = structu(i).EigenValC*PixelSize^2;

end

end