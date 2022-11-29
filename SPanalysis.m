%Author @Aiden Ross 7/8/2022

%Data input
%filename = "/Users/aidenross/Desktop/PZT_111/P-E_Loop111/40/12/ps.00000000.dat";
%filename = "PSTO_BiLayer/MartinPotential/MultiLayer/P-E_LoopV16/Polar.00034000.dat";
%filename ="/Users/aidenross/Desktop/height20_loop/height20_loop/r11h10/ps.00008000.dat";
%filename = "PSTO_BiLayer/MartinPotential/ThinFilms/P-E_LoopV7/0.50/Polar.00030000.dat";
filename ="/PSTO_BiLayer/MartinPotential/ThinFilms/P-ELoopInterp2/Polar.000710*0.dat"

fullPath = filename;
filenames = dir(fullPath);

%%
data = getData(filename);
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

title(["Orthographic Projection Pb_0_._5Sr_0_._5TiO_3"])


c = colorbar;
c.Label.String = 'log(Number Density)';
xlabel("P1")
ylabel("P2")
xlim([-0.7, 0.7])
ylim([-0.7, 0.7])


%caxis([min(gibbs,[],'all'),min(gibbs,[],'all')+2e7])


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