function [FGIT2,FGIT2_mag] = Functional_Gradient_Index_Calculator_Ver3(T2map_cart)
%Short comings - inverse mapping too many points
% edges being forced to horizontal positions - dont know how to avoid.
%% Finding first non 0 column location for whole image
cd('/Users/rlwilson77/Desktop/MRI Codes');
T2map0 = T2map_cart;
T2map0(isnan(T2map_cart)) = 0;
[colloc,rowloc] = femur_edge_detection(T2map0,1);

flip_upsidedown_T2 = flipud(T2map0);
[collocud,rowlocud] = femur_edge_detection(flip_upsidedown_T2,0); %0 might be an issue
temp_plot = zeros(384,384);
for i=1:length(collocud)
 temp_plot(rowlocud(i),collocud(i)) = 1;
end
temp_plot = flipud(temp_plot);

[collocud,rowlocud] = femur_edge_detection(temp_plot,1);
 
%clearvars temp_plot flip_upsidedown_T2

% Creating the Femur edge map

Fmap = zeros(384,384); %Preallocate the femur edge map
for i = 1:length(colloc)
    coli = colloc(i);
    rowi = rowloc(i);
    Fmap(rowi,coli) = 1;
end

clearvars rowi coli i
%% Curve Fitting

[pp,der] = Fourth_Order_BSpline(colloc,rowloc);

% This section gets rid of "knot doubles from the curve fitting. Necessary
% for the interparc code to have unique values.
[b,m1,~] = unique(pp(1,:),'first');
[~,d1] =sort(m1);
b = b(d1);
btemp = pp(2,:);
b1 = btemp(m1);
% interparc creates evenly spaced points across the bspline
t = [0:.001:1]; % dictates the "point resolution". The more points the higher the resolution.
[pp,~,~] = interparc(t,b,b1);
pp = pp';
derx = fnval(der,pp(1,:)); %the slope of the bspline at all x locations.
%% Calculating Unit Vector Normals
% Uses the fit 4th order B-spline points (pp) and the derivative slopes
% (derx) to calculate the normal line one unit vector out from the fit
% line (un) by using the point slope formula (y2-y1)=m(x2-x1)
%
nm = zeros(size(pp,1),size(pp,2)); %a "point" on the normal line from the Bspline
unit_vector = zeros(size(pp,1),size(pp,2)); % the unit vector
un = zeros(size(pp,1),size(pp,2)); %the final normal vector with position

for i=1:length(pp)
y1 = pp(2,i);
m = (derx(1,i));
x1 = pp(1,i);
x2 = pp(1,i)-1;
y2 = (-1/m)*(x2-x1)+y1;
nm(:,i) = [x2 y2];
%unit vector equation=(y2-y1)/sqrt((x2-x1)^2+(y2-y1)^2)
sqt = sqrt((pp(1,i)-nm(1,i))^2+(pp(2,i)-nm(2,i))^2);
unit_vector(:,i) = [(nm(1,i)-pp(1,i))/sqt, (nm(2,i)-pp(2,i))/sqt];
%when the slope of the line is 0 the line is either perfectly vertical 
%or horizontal and the unit vector is 0 and NaN. The following accounts for
%that. (Could also check the unit vector variable itself. Might make for a
%better check in the future).
%
if derx(i) > 0
    un(:,i) = pp(:,i)+unit_vector(:,i);
elseif derx(i) < 0
    un(:,i) = pp(:,i)-unit_vector(:,i);
elseif derx(i) == 0 && i == 1
    un(1,i) = pp(1,i)-1;
    un(2,i) = pp(2,i);
    unit_vector(1,i) = -1;
    unit_vector(2,i) = 0;
elseif derx(i) == 0 && i == length(pp)
    un(1,i) = pp(1,i)+1;
    un(2,i) = pp(2,i);
    unit_vector(1,i) = 1;
    unit_vector(2,i) = 0;
else
    un(1,i) = pp(1,i)-unit_vector(1,i);
    un(2,i) = pp(2,i)+1;
    unit_vector(1,i) = 0;
    unit_vector(2,i) = 1;
end
end
%% Calculating "Point Values"

[steps,pu] = Point_Value_CalculatorVer2(T2map0,pp,collocud,rowlocud,unit_vector,derx);
%% Calculating FGIs

[FGI]= finite_difference_scheme(steps);

final_image_information = steps(:,1:2,:);
final_image_information(:,3,:) = FGI;

clearvars FGI
%% Fitting FGI's Back to Original T2map Locations
% Need to create single columns of position and intensity data.
% Additionally, need to get rid of all NaN values.
stepscol = final_image_information(:,1,:);
stepscol = squeeze(stepscol);
stepsrow = final_image_information(:,2,:);
stepsrow = squeeze(stepsrow);
temp = [stepscol(:,1);stepscol(:,2)];
temp1 = [stepsrow(:,1);stepsrow(:,2)];

for i=3:(size(stepscol,2))
    
   stepscol1 = [temp;stepscol(:,i)];
   stepsrow1 = [temp1;stepsrow(:,i)];
   temp = stepscol1(:,1);
   temp1 = stepsrow1(:,1);
end

stepsFGI = final_image_information(:,3,:);
stepsFGI = squeeze(stepsFGI);
temp = [stepsFGI(:,1);stepsFGI(:,2)];
for i=3:(size(stepsFGI,2))
    
   stepsFGI1 = [temp;stepsFGI(:,i)];
   temp = stepsFGI1(:,1);
end

[rowt, ~] = find(isnan(stepsFGI1));
temp = stepsFGI1;
temp(rowt)= [];
stepsFGI1 = temp;
temp = stepscol1;
temp(rowt)= [];
stepscol1 = temp;
temp = stepsrow1;
temp(rowt)= [];
stepsrow1 = temp;
%%
% Fitting the data

% set fitting options
filt = 'cubicinterp';
% fitting
[FGI_fit] = fit([stepsrow1(:,1),stepscol1(:,1)], stepsFGI1(:,1), filt);
%% Matching Corresonding rows and columns to fit data
[row, col, ~] = find(T2map0);
FGIT2 = zeros(size(T2map0,1),size(T2map0,2)); %Reconstructed FGI Map but for now just T2
for i=1:size(row,1)
    FGIT2(row(i),col(i)) = FGI_fit(row(i),col(i));
end

FGIT2(isnan(FGIT2))=0;
FGIT2_mag = abs(FGIT2);
%% Plotting Figures
maxc = max(max(FGIT2));
minc = min(min(FGIT2));
figure
im2 = imagesc(FGIT2);
colormap(jet);

im2.AlphaDataMapping = 'scaled';
im2.AlphaData = [FGIT2~=0];
axis square
c = colorbar;
c.Limits = [0,15];
%% Unit Vectorsx3
% figure, hold on
% im = imagesc(T2map0);
% im.AlphaDataMapping = 'scaled';
% im.AlphaData = [T2map0~=0];
% colormap(jet);
% for j=1:4%size(steps,3)
% for i=1:size(steps,1) 
% scatter(steps(i,1,j),steps(i,2,j),10,'r','filled');
% end
% end
% set(gca, 'YDir','reverse');
% axis square
