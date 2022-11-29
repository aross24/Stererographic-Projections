%Author @Aiden Ross 10/11/2022

%1. Import the data
filename1 = "../PSTO_BiLayer/MartinPotential/ThinFilms/P-ELoopInterp2/Polar.*.dat";
filename1 = dir(filename1);
filename = "PSTO_BiLayer/MartinPotential/ThinFilms/90_Interp2/Polar.*dat";

filenames = dir(filename);

%filenames = [filename1(6:7);filenames(1:210)];


%% set up arrays for storage
filename=fullfile(filenames(60).folder, filenames(60).name);
data = getData(filename);


X = [data(:,4), data(:,5), data(:,6)];
sz = size(X);

fulldata = zeros(length(filenames), sz(1), sz(2) );
%%
xvec =  zeros(sz(1), sz(2) ); yvec =  zeros(sz(1), sz(2) );
yvec(:,2)=1;
xvec(:,1)=1;

mag = sqrt(data(:,4).^2+data(:,5).^2+data(:,6).^2);
angX = abs(dot(data(:,4:6), yvec, 2) ./ mag);
angY = abs(dot(data(:,4:6), xvec, 2) ./ mag);
buffer = 1; % 2.5° in radians 0.02
%index = (abs(1-angX) < buffer) | (abs(1-angY) < buffer);
index = (abs(1-angX) < buffer);
data = data(index,:);

X = [data(:,4), data(:,5), data(:,6)];
sz = size(X);

fulldata = zeros(length(filenames), sz(1), sz(2) );
stepData = zeros(length(filenames),1);



%%

v= VideoWriter(['SPnoFilter.mp4'],'MPEG-4')

v.FrameRate=30;
open(v)

for i = 1:length(filenames)

filename=fullfile(filenames(i).folder, filenames(i).name);

data = getData(filename);
%data = data(index,:);

newStr=split(filenames(i).name, '.');
step = str2double(newStr{2});
stepData(i)=step;


fulldata(i,:,:)=[data(:,4), data(:,5), data(:,6)];

X = [data(:,4), data(:,5)];


ang1=0;
ang2 =0;
ang3=0;

[v_R,v_T,v_O]= getRTO(ang1,ang2,ang3);

%plot(v_O(:,1).*Pscale, -v_O(:,2).*Pscale, 'rs', 'MarkerSize',10, 'MarkerFaceColor','r')
%plot(v_T(:,1).*Pscale, -v_T(:,2).*Pscale, 'go', 'MarkerSize',10, 'MarkerFaceColor','g')
%caxis([-5e7, -3e7]);

%Density Plot Processing
resolution = 100;
[values, centers]=hist3(X, 'Nbins', [resolution,resolution]);


%Plotting             

hold on
%grid on 
axis equal
ax = gca;

    ax.FontWeight = 'Bold';
    ax.FontSize = 18;
    ax.LineWidth = 1;

c = colorbar;

%PLOTS THE POINTs
plot(data(:,4), data(:,5), 'k.')

%PLOTS THE DENSITY

imagesc(centers{:}, log(values.'), 'AlphaData', log(values.'))

%imagesc(centers{:}, (values.'))


%Labels etc.

title("Orthographic Projection Pb_0_._5Sr_0_._5TiO_3")


c = colorbar;
c.Label.String = 'log(Number Density)';
caxis([0,8])
xlabel("P1")
ylabel("P2")
xlim([-0.7, 0.7])
ylim([-0.7, 0.7])
title(i)
frame = getframe(gcf);
writeVideo(v, frame)
hold off
clf
hold on

end

close(v)
%%
hold on
startStep = 60;
endStep =128;
X=[fulldata(startStep,:,1)', fulldata(startStep,:,2)'];
resolution = 100;
[values, centers]=hist3(X, 'Nbins', [resolution,resolution]);
%imagesc(centers{:}, log(values.'), 'AlphaData', log(values.'))

for i=1:11:sz(1)
    %plot(fulldata(startStep,i,1), fulldata(startStep,i,2), 'ob','MarkerSize',2)
    %plot(fulldata(endStep,i,1), fulldata(endStep,i,2), 'ok', 'MarkerSize',2)

    patch([fulldata(startStep:endStep,i,1); nan], [fulldata(startStep:endStep,i,2) ;nan], [(stepData(startStep:endStep))-startStep*1000; nan],'LineWidth',0.5,'edgecolor', 'interp', 'edgealpha',0.05);

end
%imagesc(centers{:}, log(values.'), 'AlphaData', log(values.'))

hold on
box on
axis equal
ax = gca;

    ax.FontWeight = 'Bold';
    ax.FontSize = 18;
    ax.LineWidth = 1;

xlabel("P1")
ylabel("P2")
xlim([-0.7, 0.7])
ylim([-0.7, 0.7])
colorbar;
%caxis([71350,72200])

%% Finding path density
hold on
X = zeros(sz(1)*length(filenames) , 2);

for i=1:length(filenames)

   %plot(fulldata(30:end,i,1), fulldata(30:end,i,2), '.')

   resolution = 256;
   if i >1
    lowBnd=i*(sz(1));
    upBnd =(i+1)*(sz(1))-1;
   else
      lowBnd=1;
      upBnd =(sz(1));
   end

   X(lowBnd:upBnd,:) = [fulldata(i,:,1);fulldata(i,:,2)]';




end
   [values, centers]=hist3(X, 'Nbins', [resolution,resolution]);
imagesc(centers{:}, log(values.'), 'AlphaData', log(values.'))

 function [p1c, p2c, p3c] = vecTransform(p1,p2,p3, ang1, ang2, ang3)


    T11 =  cosd(ang3) * cosd(ang1) - cosd(ang2) * sind(ang1) * sind(ang3);
    T12 =  cosd(ang3) * sind(ang1) + cosd(ang2) * cosd(ang1) * sind(ang3);
    T13 =  sind(ang3) * sind(ang2)       ;
    T21 = -sind(ang3) * cosd(ang1) - cosd(ang2) * sind(ang1) * cosd(ang3);
    T22 = -sind(ang3) * sind(ang1) + cosd(ang2) * cosd(ang1) * cosd(ang3);
    T23 =  cosd(ang3) * sind(ang2);
    T31 =  sind(ang2) * sind(ang1);
    T32 = -sind(ang2) * cosd(ang1);
    T33 =  cosd(ang2);

    p1c = T11.*p1 + T12.*p2 + T13.*p3;

    p2c = T21.*p1 + T22.*p2 + T23.*p3;

    p3c = T31.*p1 + T32.*p2 + T33.*p3;

end

function [p1l, p2l, p3l]=vecInv(p1g, p2g, p3g, ang1, ang2, ang3)


    T = zeros(3,3);
    T(1,1) =  cosd(ang3) * cosd(ang1) - cosd(ang2) * sind(ang1) * sind(ang3);
    T(1,2) =  cosd(ang3) * sind(ang1) + cosd(ang2) * cosd(ang1) * sind(ang3);
    T(1,3) =  sind(ang3) * sind(ang2)       ;
    T(2,1) = -sind(ang3) * cosd(ang1) - cosd(ang2) * sind(ang1) * cosd(ang3);
    T(2,2) = -sind(ang3) * sind(ang1) + cosd(ang2) * cosd(ang1) * cosd(ang3);
    T(2,3) =  cosd(ang3) * sind(ang2);
    T(3,1) =  sind(ang2) * sind(ang1);
    T(3,2) = -sind(ang2) * cosd(ang1);
    T(3,3) =  cosd(ang2);

    invT = T';

    p1l = invT(1,1).*p1g + invT(1,2).*p2g + invT(1,3).*p3g;

    p2l = invT(2,1).*p1g + invT(2,2).*p2g + invT(2,3).*p3g;

    p3l = invT(3,1).*p1g + invT(3,2).*p2g + invT(3,3).*p3g;


end


function [v_R,v_T,v_O]= getRTO(ang1,ang2,ang3)

RTO_directions = zeros(26,3);

RTO_directions(1,1:3)  = [+1 +1 +1]; %#1%%R1+ 
RTO_directions(2,1:3)  = [-1 -1 -1]; %#2%%R1-
RTO_directions(3,1:3)  = [-1 +1 +1]; %#3%%R2+
RTO_directions(4,1:3)  = [+1 -1 -1]; %#4%%R2-
RTO_directions(5,1:3)  = [-1 -1 +1]; %#5%%R3+
RTO_directions(6,1:3)  = [+1 +1 -1]; %#6%%R3-
RTO_directions(7,1:3)  = [+1 -1 +1]; %#7%%R4+
RTO_directions(8,1:3)  = [-1 +1 -1]; %#8%%R4-
RTO_directions(9,1:3)  = [+1 +1 0];  %#9%%O1+
RTO_directions(10,1:3) = [-1 -1 0];  %#10%%O1-
RTO_directions(11,1:3) = [+1 -1 0];  %#11%%O2+
RTO_directions(12,1:3) = [-1 +1 0];  %#12%%O2-
RTO_directions(13,1:3) = [+1 0 +1];  %#13%%O3+
RTO_directions(14,1:3) = [-1 0 -1];  %#14%%O3-
RTO_directions(15,1:3) = [+1 0 -1];  %#15%%O4+
RTO_directions(16,1:3) = [-1 0 +1];  %#16%%O4-
RTO_directions(17,1:3) = [0 +1 +1];  %#17%%O5+
RTO_directions(18,1:3) = [0 -1 -1];  %#18%%O5-
RTO_directions(19,1:3) = [0 +1 -1];  %#19%%O6+
RTO_directions(20,1:3) = [0 -1 +1];  %#20%%O6-
RTO_directions(21,1:3) = [+1 0 0];   %#21%%T1+
RTO_directions(22,1:3) = [-1 0 0];   %#22%%T1-
RTO_directions(23,1:3) = [0 +1 0];   %#23%%T2+
RTO_directions(24,1:3) = [0 -1 0];  % %#24%%T2-
RTO_directions(25,1:3) = [0 0 +1];   %#25%%T3+
RTO_directions(26,1:3) = [0 0 -1];   %#26%%T3-

v_R = RTO_directions(1:8,:)./sqrt(3);
v_O = RTO_directions(9:20,:)./sqrt(2);
v_T = RTO_directions(21:26,:);

[v_R(:,1), v_R(:,2), v_R(:,3)] = vecTransform(v_R(:,1), v_R(:,2), v_R(:,3), ang1, ang2, ang3);

[v_O(:,1), v_O(:,2), v_O(:,3)] = vecTransform(v_O(:,1), v_O(:,2), v_O(:,3), ang1, ang2, ang3);

[v_T(:,1), v_T(:,2), v_T(:,3)] = vecTransform(v_T(:,1), v_T(:,2), v_T(:,3), ang1, ang2, ang3);
v_R = v_R((v_R(:,3)>0),:);
v_O = v_O((v_O(:,3)>0),:);
v_T = v_T((v_T(:,3)>0),:);
end

function data = getData(filename)

%imports the ps.dat file into the program 
input = importdata(filename,' ',1);


header = input.textdata;
headerArray = cell2mat(header);
headerstr = convertCharsToStrings(headerArray);
C = strsplit(headerstr);
data = input.data;
nx = str2num(C(2));
ny = str2num(C(3));
nz = str2num(C(4));
data = data(((data(:,4).^2+data(:,5).^2+data(:,6).^2)>0),:);

end

function saveData(data,fileName)
% SAVEDATA saves the data into a .dat file

disp('Saving Data')
%Finding nx, ny, and nz

nx = max(data(:,1));
ny = max(data(:,2));
nz = max(data(:,3));

%creating the first line of the file
fileID = fopen(strcat(fileName,'.dat'),'w');
fprintf(fileID,'%6i %6i %6i\n',nx,ny,nz);
%writing the rest of the file
for i = 1: length(data)
    fprintf(fileID,'%6i %6i %6i %9.5f %9.5f %9.5f\n',data(i,1), data(i,2), data(i,3), data(i,4), data(i,5), data(i,6));
end
disp('File Saved')
end