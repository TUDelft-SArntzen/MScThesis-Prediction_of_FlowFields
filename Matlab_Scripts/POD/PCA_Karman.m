%%% POD-Emulator script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stefan Arntzen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% V2 - 26.02.2019 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
close all
clear all
addpath('F:\MThesis\UNIX\Documents\Thesis\5.CODE\Figures_report')
%% Readloop Tables: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic 
% Data-location 2D-Karman-street %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('F:\MThesis\UNIX\Desktop\CFD\2D_KarmanDOE\DOE\Param=2.00\DataTables')
% cd('F:\MThesis\UNIX\Desktop\CFD\2D_Karman\Case0\Karman_RoughGRID\Data_Tables\10000_12710\Data_files');
cd('F:\MThesis\UNIX\Desktop\CFD\2D_Karman\Case1\Data_Tables\8000_9999');
% Registering the amount of data-tables
dummy = dir('*.csv');
numfiles = numel(dummy);

% Initialize - preallocating
mydata = cell(1, numfiles);
Number_of_columns=10;                                       %manual entry
Number_of_gridpoints=3152;                                 %manual entry

data_total = [];
C=zeros(numfiles,Number_of_columns);

% Read loop 
for k = 1:numfiles
myfilename = dummy(k).name; 
disp(myfilename)
    
data_old = data_total;
data_new = csvread(myfilename,1,0);
data_total = horzcat(data_old, data_new);

     for k2=1:Number_of_columns
     C(k,k2) = (k2-Number_of_columns) +Number_of_columns.*k; %Finding the right index for each flow property
     end
end
Pressure   =data_total(:,C(:,2));
Velo_Mag   =data_total(:,C(:,3));
Velo_X     =data_total(:,C(:,4));
Velo_Y     =data_total(:,C(:,5));
Velo_Z     =data_total(:,C(:,6));
Time       =data_total(:,C(:,7));       % THIS CHANGES - check if correct!
% Vorticity  =data_total(:,C(:,2));

X_coor = data_total(:,C(:,8));
Y_coor = data_total(:,C(:,9));
% Z_coor = data_total(:,C(:,10));
%  clearvars -except Coefficients_c1 Coefficients_c1_5 data_total Vorticity Velo_X Velo_Y Velo_Z Time X_coor Y_coor Z_coor Pressure Temperature Time Velo_Mag Vorticity
toc
%% POD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% POD Description (brief)
%%% Modes = prinicipal component coefficients (eigenmodes/eigenvectors)
%%% Coefficients = Temporal coefficients (Rows = time, Columns = modes)
%%% Variance = Eigenvalues of the covariance matrix
%%% tsquared = sum of squares of the standarized scores for each time
%%% explained = percentage of total variance captured for each mode
%%% mu = extimated means of the variable - MEAN MODE 
tic
[Modes,Coefficients,Variance,tsquared, Explained, mu]=pca(Vorticity');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Variance - "amount of energy" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Total_variance = sum(Variance);
dummy=[0; cumsum(Explained)];                        % dummy values to add Y=0
dummy2=[0,linspace(1,size(Modes,2),size(Modes,2))];  % dummy values to add X=0

figure(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScriptAspectFigures_Pre_Report;
hold on
plot(Explained,'o-','MarkerSize',5)
plot(dummy2,dummy,'o-','MarkerSize',5)
legend('Variance','Cumulative variance','Location', 'Northeast')
xlim([0 60]); 
% title('Pressure variance plot');
xlabel('Number of modes');ylabel('Variance [%]')
% parameters for saving the figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figName = 'Variance_vorticity';
% figName = 'variance_pressure';
Loop=false;
saveFig =true;    %true = safe, false = don't safe
Extension = true;   %true = PDF , false = PNG
ScriptAspectFigures_Post_Report;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data reduction - Mode truncation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Modes = Modes(:,1:50);
Coefficients = Coefficients(:,1:50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MESH generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=X_coor(:,1); 
Y=Y_coor(:,1);
% dTime=Time(1,2)-Time(1,1);                                 %delta Time
tri=delaunay(X,Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Deterimine position probe points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPoint_suitable = find( abs(X-0.048)<5.1e-4 & abs(Y)<5.1e-4);
PP001=find(abs(X-0.0481)<5.1e-4 & abs(Y)<5.1e-4);      Pressure_PP001=Pressure(PP001(1),:);
PP050=find(abs(X-0.055)<5e-3 & abs(Y)<5.1e-4);         Pressure_PP050=Pressure(PP050(1),:);
PP200=find(abs(X-0.06)<5e-3 & abs(Y)<5.1e-4);          Pressure_PP200=Pressure(PP200(1),:);
PP650=find(abs(X-0.075)<5e-3 & abs(Y)<5.1e-4);         Pressure_PP650=Pressure(PP650(1),:);
PP950=find(abs(X-0.095)<5e-3 & abs(Y)<5.1e-4);         Pressure_PP950=Pressure(PP950(1),:);

% figure; hold on
% triplot(tri,X,Y)
% plot(X(PP001(1)),Y(PP001(1)),'r.','MarkerSize',20)
% plot(X(PP050(1)),Y(PP050(1)),'r.','MarkerSize',20)
% plot(X(PP200(1)),Y(PP200(1)),'r.','MarkerSize',20)
% plot(X(PP650(1)),Y(PP650(1)),'r.','MarkerSize',20)
% plot(X(PP950(1)),Y(PP950(1)),'r.','MarkerSize',20)
% axis equal tight

% Clear unwanted variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dummy dummy2
%% Flow-field simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocating the matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reconstructed_mode=zeros(size(Modes,1),size(Modes,2))';
Flow_field=zeros(1,size(Modes,1));
Fluc_flow_field=zeros(1,size(Modes,1));
% Pressure_point_1801_POD=zeros(1,size(Time,2));
% Pressure_cl_pod=zeros(size(Centerline,1),size(Time,2));
%%% Mode 0, mean mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mean_Field=mu';
% figure
% trisurf(tri,X,Y,Mean_Field)
% % viscircles([0.04 0],0.00125,'Color','k')
% colormap(jet), hcb=colorbar;
% % caxis([20 70]) %velocity colour axis
% % caxis([-10150 -9500])
% % caxis([0 2000])
% shading interp, view(2), axis equal tight
% xlim([0 1.2]);ylim([-0.1 0.1])
% ylabel(hcb, 'Vorticity (/s)')
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 18 5.5]);
% xlabel('X (m)');ylabel('Y (m)');title('Mean mode - Case 1')
%%%% Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for temploop=1:size(Time,2)
Time_sim = Time(1,1) + temploop*dTime;                            %simulation time 

    for m=1:size(Modes,2)
        Reconstructed_mode(m,:) = Coefficients(temploop,m)' * Modes(:,m);
    end    
  
%%% Mean mode + Fluctuations 
Flow_field(temploop,:) = Mean_Field' + sum(Reconstructed_mode,1);
Fluc_flow_field(temploop,:) = sum(Reconstructed_mode,1);

%%% Niet nodig ---> probe point van POD kan gewoon buiten de for loop.
%%% Pressure point -> fft  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Find point which is suited for FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PPoint_suitable = find( abs(X-0.048)<5.1e-4 & abs(Y)<5.1e-4);
% Pressure_point_1801 = Pressure(1801,:);
% Pressure_point_1801_POD(1,t) = Flow_field(1,1801); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Centerline data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure_cl_pod(:,t) = Flow_field(1, Centerline) - mean(Flow_field(1, Centerline) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Flow-field figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% viscircles([0.04 0],0.00125,'Color','k');
% caxis([-10500 -9500]); 
% % 
% trisurf(tri,X,Y,Fluc_flow_field(t,:)');
% colormap(jet), hcb=colorbar;
% shading interp, view(2), axis equal tight
% xlabel('X (m)');ylabel('Y (m)')
% xlim([0 0.1]); caxis([-200 200])
% viscircles([0.04 0],0.00125,'Color','k');
% title(sprintf('Mode 1 to 50, (Time = %d s)',Time_sim))
% drawnow
% ylabel(hcb, 'Absolute pressure (Pa)')
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 18 5.5]);
% print('-dpng','-r0',sprintf('image%d.png',t));
end
%% Interpolation temporal Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Attempt 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Time domain interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Coefficients_ip = 0.5*(Coefficients_c1(:,:) + Coefficients_c0(1:1000,:));

Reconstructed_mode_ip=zeros(size(Modes_c1,2),size(Modes_c1,1));
Flow_field_ip=zeros(size(Modes_c1,2),size(Modes_c1,1));
% X=X_coor(:,1); 
% Y=Y_coor(:,1);
% tri=delaunay(X,Y);
mm=1000;
%%%% Time loop
for temploop=206
  for m=1:size(Modes_c1,2)
        Reconstructed_mode_ip(m,:) = Coefficients_ip(temploop,m) * Modes_c0(:,m)';
  end
  
Flow_field_ip(temploop,:) =  sum(Reconstructed_mode_ip);

Loop=true;
ScriptAspectFigures_Pre_Report;
trisurf(tri_c0,X_c0*mm,Y_c0*mm,Flow_field_ip(temploop,:)');
viscircles([0.04*mm 0],0.00125*mm,'Color','k');
colormap(jet), hcb=colorbar;
shading interp, view(2), axis equal tight
xlabel('X (mm)');ylabel('Y (mm)')
xlim([0.03*mm 0.1*mm]); ylim([-0.01*mm 0.01*mm]);
caxis([-200 200])
drawnow
ylabel(hcb, 'Pressure (Pa)')
ScriptAspectFigures_Post_Report;

% print('-dpng','-r0',sprintf('image%d.png',t));
  
% subplot(2,1,1)
% plot(Time(1,[1:1000]),Flow_field_c01(:,1),'r-*')
% legend('Mode 1')
% ylim([-2300 2300])
% title('Case interpolated, Mass flow 1.325g/s');xlabel('Time [s]');ylabel('POD-coefficients [-]')
% grid on, grid minor
% subplot(2,1,2)
% trisurf(tri,X,Y,Flow_field_c01(t,:)');
% viscircles([0.04 0],0.00125,'Color','k')
% caxis([-70  70]), 
% colormap(jet), colorbar
% shading interp, view(2), axis equal tight
% xlabel('Time [s]');ylabel('Y-coordinate [m]')
% % % title(sprintf('Mode 1, interpolated case, (Time = %d s)',Time_sim))
% drawnow
% print('-dpng','-r0',sprintf('image%d.png',t));
% 
end
PPOD_001_ip = Flow_field_ip(:,PP001(1)); 
PPOD_050_ip = Flow_field_ip(:,PP050(1));
PPOD_200_ip = Flow_field_ip(:,PP200); 
PPOD_650_ip = Flow_field_ip(:,PP650);
PPOD_950_ip = Flow_field_ip(:,PP950); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Attempt 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Frequency domain interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Coefficients=Coefficients([1:999],[1:998]);
% Coefficients_fft=fft(Coefficients(:,3));                %Score Case0- FFT
% Coefficients_fft_C11=fft(Coefficients_C11_temp(:,3));   %Score Case1- FFT     
% % score_fft_c0_temp=score_fft_c0([1:999],[1:998]); %Data-truncation - same size
% % score_fft_c01=(score_fft_c0_temp + score_fft_c1)/2;    % Interpolation in freq domain
% % 
% % score_c01=ifft(score_fft_c01);                   %Inverse FFT back to time-domain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Attempt 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Artificial signal recreation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temp = 1/Fs;
% temp = (0:Lts-1)*Temp;

% S1 = 1891*sin(2*pi*1870*temp);
% S2 = 593.6*sin(2*pi*1870*temp);
% S3 = 187.5*sin(2*pi*3600*temp);
% Y1= fft(S1);
% 
% P2temp = abs(Y1/Lts);
% P1temp = P2temp(1:Lts/2+1);
% P1temp(2:end-1) = 2*P1temp(2:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Pressure spectra for recreation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure
% % plot(f,P1,f,P1temp,f,P3)
% % title('Single-Sided Amplitude Spectrum of Beta(t)')
% % xlabel('f (Hz)')
% % ylabel('|Magnitude(f)|')
% % grid on, grid minor
% % xlim([0 5000])
% % legend('Case 1','Case ip','Case 1.1')
%%%% Temporal recreated signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure
% % plot(temp,S1,'bo-',temp,S2,'r*-',temp,S3,'gx-')
% % legend('Case ip, mode 1- beta signal','Case ip, mode 2- beta signal','Case ip, mode 3- beta signal')
% % xlabel('X-coor [m]');ylabel('Beta-coeff');
% % grid on, grid minor
%%%% Mesh generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X=X_coor(:,1); 
% Y=Y_coor(:,1);
% tri=delaunay(X,Y);
% tri2=delaunay(X_C11_temp,Y_C11_temp); 

% temptemp=Modes(:,1);
% temptemp2=Modes_C11_temp(:,1);
% temptemp3=1/2*(temptemp+temptemp2);
% 
% temptemp_2=Modes(:,2);
% temptemp_2_2=Modes_C11_temp(:,2);
% temptemp_2_3=1/2*(temptemp_2+temptemp_2_2);
% 
% temptemp_3=Modes(:,3);
% temptemp_3_2=Modes_C11_temp(:,3);
% temptemp_3_3=1/2*(temptemp_2+temptemp_2_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flow field reconstruction for inerpolated modes + recreated signals %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for temptemptemptemp=1:size(S1,2)
% flowfield_temp=mu'+(S1(1,temptemptemptemp)*temptemp3) +(S2(1,temptemptemptemp)*temptemp_2_3) + (S3(1,temptemptemptemp)*temptemp_3_3);
% trisurf(tri,X,Y,flowfield_temp');
% viscircles([0.04 0],0.00125,'Color','k')
% caxis([-10500 -9500]), 
% colormap(jet), colorbar
% shading interp,view(2), axis equal tight 
% xlabel('Time [s]');ylabel('Y-coordinate [m]')
% title('Mean mode Case1, Mode 1 to 3 recreated')
% % drawnow
% print('-dpng','-r0',sprintf('image%d.png',temptemptemptemp));
% end
%% Pressure & FFT plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Determine the %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPOD_001 = Fluc_flow_field(:,PP001(1)); 
PPOD_050 = Fluc_flow_field(:,PP050(1));
PPOD_200 = Fluc_flow_field(:,PP200(1)); 
PPOD_650 = Fluc_flow_field(:,PP650(1));
PPOD_950 = Fluc_flow_field(:,PP950(1)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Temporal figure of static pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure_PP001=PPOD_001_ip;
% Pressure_PP050=PPOD_050_ip;
% Pressure_PP200=PPOD_200_ip;
% Pressure_PP650=PPOD_650_ip;
% Pressure_PP950=PPOD_950_ip;

%%%
Mean_PP001=zeros(1,size(Time,2));Mean_PP050=zeros(1,size(Time,2));
Mean_PP200=zeros(1,size(Time,2));Mean_PP650=zeros(1,size(Time,2));
Mean_PP950=zeros(1,size(Time,2));

for d=1:size(Time,2)
    Mean_PP001(d)=mean(Pressure_PP001(1:d));
    Mean_PP050(d)=mean(Pressure_PP050(1:d)); 
    Mean_PP200(d)=mean(Pressure_PP200(1:d)); 
    Mean_PP650(d)=mean(Pressure_PP650(1:d)); 
    Mean_PP950(d)=mean(Pressure_PP950(1:d)); 
end
%%% temporal figure of pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lwdth=1.25;
saveFig = false;
Extension =false;   %true = PDF , false = PNG

ScriptAspectFigures_Pre;
hold on
plot(Time(1,:),Pressure_PP001,'-','LineWidth',Lwdth)
% plot(Time(1,:),Pressure_PP050,'-','LineWidth',Lwdth)
plot(Time(1,:),Pressure_PP200,'-','Color',[0.9290, 0.6940, 0.1250],'LineWidth',Lwdth)
% plot(Time(1,:),Pressure_PP650,'-','LineWidth',Lwdth)
% plot(Time(1,:),Pressure_PP950,'-','LineWidth',Lwdth)

ax = gca;
ax.ColorOrderIndex = 1;

 plot(Time(1,:),Mean_PP001,'--','LineWidth',Lwdth)
% plot(Time(1,:),Mean_PP050,'--','LineWidth',Lwdth)
plot(Time(1,:),Mean_PP200,'--','Color',[0.9290, 0.6940, 0.1250],'LineWidth',Lwdth)
% plot(Time(1,:),Mean_PP650,'--','LineWidth',Lwdth)
% plot(Time(1,:),Mean_PP950,'--','LineWidth',Lwdth)
%%% axis properties and labels 
temptitle=('Pressure measured at probe points');
xlabel('Time (s)'),ylabel('Pressure (Pa)');
xlim([Time(1,1)-dTime Time(1,end)+dTime]),
% ylim([-200 250])

% legend('Probe point @ X=1','Probe point @ X=50','Probe point @ X=200','Probe point @ X=650'...
%     ,'Probe point @ X=950','Location', 'best');
legend('P.P @ X=1 mm','P.P @ X=200 mm','Mean P.P @ X=1 mm','Mean P.P @ X=200 mm','Location', 'northeast')
set(legend,'NumColumns',2)
hold off
figName = 'TemporalPressure_meanPressure';
ScriptAspectFigures_Post;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% LES vs Truncated POD figures..%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(Time(1,:),Pressure_PP001,'-','LineWidth',Lwdth)
plot(Time(1,:),PPOD_001,'--','LineWidth',Lwdth)
xlabel('Time (s)'),ylabel('Pressure (Pa)')
xlim([Time(1,1)-dTime Time(1,end)+dTime])
% ylim([-200 250])
legend('LES-data','Truncated POD-data','Location','northeast')
set(legend,'NumColumns',2)
temptitle=('Pressure signal for P.P at X=1 mm');
figName = 'PressureSignal_X1_LESvsPOD';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(Time(1,:),Pressure_PP050,'-','LineWidth',Lwdth)
plot(Time(1,:),PPOD_050,'--','LineWidth',Lwdth)
xlabel('Time (s)'),ylabel('Pressure (Pa)')...
    ,xlim([Time(1,1)-dTime Time(1,end)+dTime]);
% ylim([-200 250])
legend('LES-data','Truncated POD-data','Location','northeast')
temptitle=('Pressure signal for P.P at X=50 mm');
set(legend,'NumColumns',2)
figName = 'PressureSignal_X50_LESvsPOD';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(Time(1,:),Pressure_PP200,'-','LineWidth',Lwdth)
plot(Time(1,:),PPOD_200,'--','LineWidth',Lwdth)
xlabel('Time (s)'),ylabel('Pressure (Pa)')...
    ,xlim([Time(1,1)-dTime Time(1,end)+dTime]);
% ylim([-200 250])
legend('LES-data','Truncated POD-data','Location','northeast')
temptitle=('Pressure signal for P.P at X=200 mm');
set(legend,'NumColumns',2)
figName = 'PressureSignal_X200_LESvsPOD';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(Time(1,:),Pressure_PP650,'-','LineWidth',Lwdth)
plot(Time(1,:),PPOD_650,'--','LineWidth',Lwdth)
xlabel('Time (s)'),ylabel('Pressure (Pa)')...
    ,xlim([Time(1,1)-dTime Time(1,end)+dTime]);
% ylim([-200 250])
legend('LES-data','Truncated POD-data','Location','northeast')
temptitle=('Pressure signal for P.P at X=650 mm');
set(legend,'NumColumns',2)
figName = 'PressureSignal_X650_LESvsPOD';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(Time(1,:),Pressure_PP950,'-','LineWidth',Lwdth)
plot(Time(1,:),PPOD_950,'--','LineWidth',Lwdth)
xlabel('Time (s)'),ylabel('Pressure (Pa)')...
    ,xlim([Time(1,1)-dTime Time(1,end)+dTime]);
% ylim([-200 250])
legend('LES-data','Truncated POD-data','Location','northeast')
set(legend,'NumColumns',2)
temptitle=('Pressure signal for P.P at X=950 mm');
figName = 'PressureSignal_X950_LESvsPOD';
ScriptAspectFigures_Post;

clear saveFig Extension Lwdth
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sampling and time specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dTime=Time(1,2)-Time(1,1);
Fs=1/dTime;                             % sampling freq
Lts = size(Time,2);                     % length time signal
f = Fs*(0:(Lts/2))/Lts;                 % freq range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DC offset subtraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dummy_DCsubtract = 1:size(Time,2)
    Pressure_PP001_DCsubtract=Pressure_PP001-mean(Pressure_PP001);
    Pressure_PP050_DCsubtract=Pressure_PP050-mean(Pressure_PP050);
    Pressure_PP200_DCsubtract=Pressure_PP200-mean(Pressure_PP200);
    Pressure_PP650_DCsubtract=Pressure_PP650-mean(Pressure_PP650);
    Pressure_PP950_DCsubtract=Pressure_PP950-mean(Pressure_PP950); 
%temp = Pressure_point_1801 - mean(Pressure_point_1801);
%temp_pod = Pressure_point_1801_POD - mean(Pressure_point_1801_POD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Standard FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PP001_fft=fft(Pressure_PP001_DCsubtract); PPOD_001_fft=fft(PPOD_001);
PP050_fft=fft(Pressure_PP050_DCsubtract); PPOD_050_fft=fft(PPOD_050);
PP200_fft=fft(Pressure_PP200_DCsubtract); PPOD_200_fft=fft(PPOD_200);
PP650_fft=fft(Pressure_PP650_DCsubtract); PPOD_650_fft=fft(PPOD_650);
PP950_fft=fft(Pressure_PP950_DCsubtract); PPOD_950_fft=fft(PPOD_950);
% P_fft_pod=fft(temp_pod);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fourier window specification + windown shift %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pds_PP001 = abs(PP001_fft/Lts);                     % double-sided amplitude spectrum
Pss_PP001 = Pds_PP001(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PP001(2:end-1) = 2*Pss_PP001(2:end-1);          % frequency shift

Pds_PP050 = abs(PP050_fft/Lts);                     % double-sided amplitude spectrum
Pss_PP050 = Pds_PP050(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PP050(2:end-1) = 2*Pss_PP050(2:end-1);          % frequency shift

Pds_PP200 = abs(PP200_fft/Lts);                     % double-sided amplitude spectrum
Pss_PP200 = Pds_PP200(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PP200(2:end-1) = 2*Pss_PP200(2:end-1);          % frequency shift

Pds_PP650 = abs(PP650_fft/Lts);                     % double-sided amplitude spectrum
Pss_PP650 = Pds_PP650(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PP650(2:end-1) = 2*Pss_PP650(2:end-1);          % frequency shift

Pds_PP950 = abs(PP950_fft/Lts);                     % double-sided amplitude spectrum
Pss_PP950 = Pds_PP950(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PP950(2:end-1) = 2*Pss_PP950(2:end-1);          % frequency shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pds_PPOD_001 = abs(PPOD_001_fft/Lts);                     % double-sided amplitude spectrum
Pss_PPOD_001 = Pds_PPOD_001(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PPOD_001(2:end-1) = 2*Pss_PPOD_001(2:end-1);          % frequency shift

Pds_PPOD_050 = abs(PPOD_050_fft/Lts);                     % double-sided amplitude spectrum
Pss_PPOD_050 = Pds_PPOD_050(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PPOD_050(2:end-1) = 2*Pss_PPOD_050(2:end-1);          % frequency shift

Pds_PPOD_200 = abs(PPOD_200_fft/Lts);                     % double-sided amplitude spectrum
Pss_PPOD_200 = Pds_PPOD_200(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PPOD_200(2:end-1) = 2*Pss_PPOD_200(2:end-1);          % frequency shift

Pds_PPOD_650 = abs(PPOD_650_fft/Lts);                     % double-sided amplitude spectrum
Pss_PPOD_650 = Pds_PPOD_650(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PPOD_650(2:end-1) = 2*Pss_PPOD_650(2:end-1);          % frequency shift

Pds_PPOD_950 = abs(PPOD_950_fft/Lts);                     % double-sided amplitude spectrum
Pss_PPOD_950 = Pds_PPOD_950(1:Lts/2+1);                   % single-sided amplitude spectrum
Pss_PPOD_950(2:end-1) = 2*Pss_PPOD_950(2:end-1);          % frequency shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Frequency domain plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFig   = true;
Extension = false;   %true = PDF , false = PNG
Lwdth             = 1.25;  %linewidth

hold on
ScriptAspectFigures_Pre;
plot(f,Pss_PP001,'-','LineWidth',Lwdth)
plot(f,Pss_PP050,'-','LineWidth',Lwdth)
plot(f,Pss_PP200,'-','LineWidth',Lwdth)
plot(f,Pss_PP650,'-','LineWidth',Lwdth)
plot(f,Pss_PP950,'-','LineWidth',Lwdth)
ylabel('Pressure (Pa)')
xlabel('Frequency (Hz)')
xlim([0 1000])
temptitle=('Pressure-spectra');
legend('P.P @ X=1','P.P @ X=50','P.P @ X=200','P.P @ X=650'...
    ,'P.P @ X=950','Location', 'northeast')
set(legend,'NumColumns',2)
figName = 'PressureSpectra';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(f,Pss_PP001,'-','LineWidth',Lwdth)
plot(f,Pss_PPOD_001,'--','LineWidth',Lwdth)
xlabel('Frequency (Hz)'),ylabel('Pressure (Pa)')...
    ,
% xlim([0 1000]), 
% ylim([0 120])
legend('LES-data','Truncated POD-data','Location','northeast')
temptitle=('Pressure spectra for P.P at X=1 mm');
hold off
figName = 'PressureSpectra_X1_LESvsPOD';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(f,Pss_PP050,'-','LineWidth',Lwdth)
plot(f,Pss_PPOD_050,'--','LineWidth',Lwdth)
xlabel('Frequency (Hz)'),ylabel('Pressure (Pa)')...
    ,
% xlim([0 1000]),
% ylim([0 120])
legend('LES-data','Truncated POD-data','Location','northeast')
temptitle=('Pressure spectra for P.P at X=50 mm');
hold off
figName = 'PressureSpectra_X50_LESvsPOD';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(f,Pss_PP200,'-','LineWidth',Lwdth)
plot(f,Pss_PPOD_200,'--','LineWidth',Lwdth)
xlabel('Frequency (Hz)'),ylabel('Pressure (Pa)')...
    ,
% xlim([0 1000]),ylim([0 120])
legend('LES-data','Truncated POD-data','Location','northeast')
temptitle=('Pressure spectra for P.P at X=200 mm');
hold off
figName = 'PressureSpectra_X200_LESvsPOD';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(f,Pss_PP650,'-','LineWidth',Lwdth)
plot(f,Pss_PPOD_650,'--','LineWidth',Lwdth)
xlabel('Frequency (Hz)'),ylabel('Pressure (Pa)')...
    ,
% xlim([0 1000]),ylim([0 120])
legend('LES-data','Truncated POD-data','Location','northeast')
temptitle=('Pressure spectra for P.P at X=650 mm');
hold off
figName = 'PressureSpectra_X650_LESvsPOD';
ScriptAspectFigures_Post;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ScriptAspectFigures_Pre;
hold on
plot(f,Pss_PP950,'-','LineWidth',Lwdth)
plot(f,Pss_PPOD_950,'--','LineWidth',Lwdth)
xlabel('Frequency (Hz)'),ylabel('Pressure (Pa)')...
    ,
% xlim([0 1000]),ylim([0 120])
legend('LES-data','Truncated POD-data','Location','northeast')
temptitle=('Pressure spectra for P.P at X=950 mm');
hold off
figName = 'PressureSpectra_X950_LESvsPOD';
ScriptAspectFigures_Post;

close all
clear Extension saveFig Lwdth
%% Error quantification%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error plots - comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dTime=Time(1,2)-Time(1,1);
tallocation=size(Time,2);
ngrid=size(X_coor,1);

%%%% define wake area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=X_coor(:,1);
Y=Y_coor(:,1);
LLimit=0.035; ULimit=0.05;

%%% find indexes for certain domains -- VON KARMAN %%%%%%%%%%%%%%%%%%%%%%%%
idx_wake=find(X>ULimit);
idx_cyl=find(X>LLimit & X<ULimit & Y<0.005 & Y>-0.005);
idx_inlet=find(X<LLimit);

%%%% We are interested in the MRE of the fluctuating field %%%%%%%%%%%%%%%%
% Pressure_fluc=Pressure-mean(Pressure,2);
Pressure_fluc=Velo_Mag-mean(Velo_Mag,2);

%%%% Allocating the matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% mode 1 to 1
Reconstructed_mode_1to1=zeros(1,ngrid);
Fluc_flow_field_m1to1=zeros(tallocation,ngrid);
Full_flow_field_m1to1=zeros(tallocation,ngrid);
% mode 1 to 5
Reconstructed_mode_1to5=zeros(5,ngrid);
Fluc_flow_field_m1to5=zeros(tallocation,ngrid);
% mode 1 to 10
Reconstructed_mode_1to10=zeros(10,ngrid);
Fluc_flow_field_m1to10=zeros(tallocation,ngrid);
% mode 1 to 50
Reconstructed_mode_1to50=zeros(50,ngrid);
Fluc_flow_field_m1to50=zeros(tallocation,ngrid);
% mode 1 to 100
Reconstructed_mode_1to100=zeros(100,ngrid);
Fluc_flow_field_m1to100=zeros(tallocation,ngrid);
% mode 1 to 250
Reconstructed_mode_1to250=zeros(250,ngrid);
Fluc_flow_field_m1to250=zeros(tallocation,ngrid);
% mode 1 to 500
Reconstructed_mode_1to500=zeros(250,ngrid);
Fluc_flow_field_m1to500=zeros(tallocation,ngrid);

tic
for t2=1:tallocation
    for m=1
        %%% Mode 1 to 1
        Reconstructed_mode_1to1(m,:) = Coefficients(t2,m)' * Modes(:,m);
        Fluc_flow_field_m1to1(t2,:) = sum(Reconstructed_mode_1to1,1);
        
%         Full_flow_field_m1to1(t2,:) = Mean_Field'+sum(Reconstructed_mode_1to1,1);
    end  
    for m2=1:5
        %%% Mode 1 to 5        
        Reconstructed_mode_1to5(m2,:) = Coefficients(t2,m2)' * Modes(:,m2);
        Fluc_flow_field_m1to5(t2,:) = sum(Reconstructed_mode_1to5,1);
    end
    for m3=1:10
        %%% Mode 1 to 10   
        Reconstructed_mode_1to10(m3,:) = Coefficients(t2,m3)' * Modes(:,m3);
        Fluc_flow_field_m1to10(t2,:) = sum(Reconstructed_mode_1to10,1);
    end
    for m4=1:50
        %%% Mode 1 to 50    
        Reconstructed_mode_1to50(m4,:) = Coefficients(t2,m4)' * Modes(:,m4);
        Fluc_flow_field_m1to50(t2,:) = sum(Reconstructed_mode_1to50,1); 
    end
    for m5=1:100
        %%% Mode 1 to 100    
        Reconstructed_mode_1to100(m5,:) = Coefficients(t2,m5)' * Modes(:,m5);
        Fluc_flow_field_m1to100(t2,:) = sum(Reconstructed_mode_1to100,1); 
    end
    for m6=1:250
        %%% Mode 1 to 250    
        Reconstructed_mode_1to250(m6,:) = Coefficients(t2,m6)' * Modes(:,m6);
        Fluc_flow_field_m1to250(t2,:) = sum(Reconstructed_mode_1to250,1); 
    end
     for m7=1:500
        %%% Mode 1 to 500    
        Reconstructed_mode_1to500(m7,:) = Coefficients(t2,m7)' * Modes(:,m7);
        Fluc_flow_field_m1to500(t2,:) = sum(Reconstructed_mode_1to500,1); 
    end
end 

%%% Domain Inlet (x<0.035)
       temp_1to1_inlet = abs(Pressure_fluc(idx_inlet,:)-Fluc_flow_field_m1to1(:,idx_inlet)');
       dev_1to1_inlet = sum(temp_1to1_inlet,1);
       % mode 1 to 5
       temp_1to5_inlet = abs(Pressure_fluc(idx_inlet,:)-Fluc_flow_field_m1to5(:,idx_inlet)');
       dev_1to5_inlet = sum(temp_1to5_inlet,1);
       % mode 1 to 10
       temp_1to10_inlet = abs(Pressure_fluc(idx_inlet,:)-Fluc_flow_field_m1to10(:,idx_inlet)');
       dev_1to10_inlet = sum(temp_1to10_inlet,1);
       % mode 1 to 50
       temp_1to50_inlet = abs(Pressure_fluc(idx_inlet,:)-Fluc_flow_field_m1to50(:,idx_inlet)');
       dev_1to50_inlet = sum(temp_1to50_inlet,1);
       % mode 1 to 100
       temp_1to100_inlet = abs(Pressure_fluc(idx_inlet,:)-Fluc_flow_field_m1to100(:,idx_inlet)');
       dev_1to100_inlet = sum(temp_1to100_inlet,1);
       % mode 1 to 250
       temp_1to250_inlet = abs(Pressure_fluc(idx_inlet,:)-Fluc_flow_field_m1to250(:,idx_inlet)');
       dev_1to250_inlet = sum(temp_1to250_inlet,1);
       % mode 1 to 500
       temp_1to500_inlet = abs(Pressure_fluc(idx_inlet,:)-Fluc_flow_field_m1to500(:,idx_inlet)');
       dev_1to500_inlet = sum(temp_1to500_inlet,1);
       
%%% Domain Cylinder (X>LLimit & X<ULimit & Y<0.005 & Y>-0.005)
       temp_1to1_cyl = abs(Pressure_fluc(idx_cyl,:)-Fluc_flow_field_m1to1(:,idx_cyl)');
       dev_1to1_cyl = sum(temp_1to1_cyl,1);
       
%        temp_fullflow=abs(Pressure(idx_cyl,:)-Full_flow_field_m1to1(:,idx_cyl)');
%        dev_fullflow =sum(temp_fullflow,1);
       % mode 1 to 5
       temp_1to5_cyl = abs(Pressure_fluc(idx_cyl,:)-Fluc_flow_field_m1to5(:,idx_cyl)');
       dev_1to5_cyl = sum(temp_1to5_cyl,1);
       % mode 1 to 10
       temp_1to10_cyl = abs(Pressure_fluc(idx_cyl,:)-Fluc_flow_field_m1to10(:,idx_cyl)');
       dev_1to10_cyl = sum(temp_1to10_cyl,1);
       % mode 1 to 50
       temp_1to50_cyl = abs(Pressure_fluc(idx_cyl,:)-Fluc_flow_field_m1to50(:,idx_cyl)');
       dev_1to50_cyl = sum(temp_1to50_cyl,1);
       % mode 1 to 100
       temp_1to100_cyl = abs(Pressure_fluc(idx_cyl,:)-Fluc_flow_field_m1to100(:,idx_cyl)');
       dev_1to100_cyl = sum(temp_1to100_cyl,1);
       % mode 1 to 250
       temp_1to250_cyl = abs(Pressure_fluc(idx_cyl,:)-Fluc_flow_field_m1to250(:,idx_cyl)');
       dev_1to250_cyl = sum(temp_1to250_cyl,1);
       % mode 1 to 500
       temp_1to500_cyl = abs(Pressure_fluc(idx_cyl,:)-Fluc_flow_field_m1to500(:,idx_cyl)');
       dev_1to500_cyl = sum(temp_1to500_cyl,1);

% %%% Domain Wake (x>0.05)
       temp_1to1_wake = abs(Pressure_fluc(idx_wake,:)-Fluc_flow_field_m1to1(:,idx_wake)');
       dev_1to1_wake = sum(temp_1to1_wake,1);
       % mode 1 to 5
       temp_1to5_wake = abs(Pressure_fluc(idx_wake,:)-Fluc_flow_field_m1to5(:,idx_wake)');
       dev_1to5_wake = sum(temp_1to5_wake,1);
       % mode 1 to 10
       temp_1to10_wake = abs(Pressure_fluc(idx_wake,:)-Fluc_flow_field_m1to10(:,idx_wake)');
       dev_1to10_wake = sum(temp_1to10_wake,1);
       % mode 1 to 50
       temp_1to50_wake = abs(Pressure_fluc(idx_wake,:)-Fluc_flow_field_m1to50(:,idx_wake)');
       dev_1to50_wake = sum(temp_1to50_wake,1);
       % mode 1 to 100
       temp_1to100_wake = abs(Pressure_fluc(idx_wake,:)-Fluc_flow_field_m1to100(:,idx_wake)');
       dev_1to100_wake = sum(temp_1to100_wake,1);
       % mode 1 to 250
       temp_1to250_wake = abs(Pressure_fluc(idx_wake,:)-Fluc_flow_field_m1to250(:,idx_wake)');
       dev_1to250_wake = sum(temp_1to250_wake,1);
       % mode 1 to 500
       temp_1to500_wake = abs(Pressure_fluc(idx_wake,:)-Fluc_flow_field_m1to500(:,idx_wake)');
       dev_1to500_wake = sum(temp_1to500_wake,1);
toc
       
%%% Base for the three domains
base_inlet=sum(abs(Pressure_fluc(idx_inlet,:)));
base_cyl=  sum(abs(Pressure_fluc(idx_cyl,:)));
base_wake= sum(abs(Pressure_fluc(idx_wake,:)));
% 
% temp_base_cyl=sum(abs(Pressure(idx_cyl,:)));
%%% allocation of MRE values
MRE_1to1_wake=zeros(1,tallocation);
MRE_1to5_wake=zeros(1,tallocation);
MRE_1to10_wake=zeros(1,tallocation);
MRE_1to50_wake=zeros(1,tallocation);
MRE_1to100_wake=zeros(1,tallocation);
MRE_1to250_wake=zeros(1,tallocation);
MRE_1to500_wake=zeros(1,tallocation);

MRE_1to1_inlet=zeros(1,tallocation);
MRE_1to5_inlet=zeros(1,tallocation);
MRE_1to10_inlet=zeros(1,tallocation);
MRE_1to50_inlet=zeros(1,tallocation);
MRE_1to100_inlet=zeros(1,tallocation);
MRE_1to250_inlet=zeros(1,tallocation);
MRE_1to500_inlet=zeros(1,tallocation);

MRE_1to1_cyl=zeros(1,tallocation);
MRE_1to5_cyl=zeros(1,tallocation);
MRE_1to10_cyl=zeros(1,tallocation);
MRE_1to50_cyl=zeros(1,tallocation);
MRE_1to100_cyl=zeros(1,tallocation);
MRE_1to250_cyl=zeros(1,tallocation);
MRE_1to500_cyl=zeros(1,tallocation);

for dummy=1:tallocation
%     temp_MRE_fullflow(:,dummy)=dev_fullflow(:,dummy)/temp_base_cyl(:,dummy);
    %%% Inlet
MRE_1to1_inlet(:,dummy)  = dev_1to1_inlet(:,dummy)/base_inlet(:,dummy);
MRE_1to5_inlet(:,dummy)  = dev_1to5_inlet(:,dummy)/base_inlet(:,dummy);
MRE_1to10_inlet(:,dummy) = dev_1to10_inlet(:,dummy)/base_inlet(:,dummy);
MRE_1to50_inlet(:,dummy) = dev_1to50_inlet(:,dummy)/base_inlet(:,dummy);
MRE_1to100_inlet(:,dummy)= dev_1to100_inlet(:,dummy)/base_inlet(:,dummy);
MRE_1to250_inlet(:,dummy)= dev_1to250_inlet(:,dummy)/base_inlet(:,dummy);
MRE_1to500_inlet(:,dummy)= dev_1to500_inlet(:,dummy)/base_inlet(:,dummy);
    %%% Cyl
MRE_1to1_cyl(:,dummy)  = dev_1to1_cyl(:,dummy)/base_cyl(:,dummy);
MRE_1to5_cyl(:,dummy)  = dev_1to5_cyl(:,dummy)/base_cyl(:,dummy);
MRE_1to10_cyl(:,dummy) = dev_1to10_cyl(:,dummy)/base_cyl(:,dummy);
MRE_1to50_cyl(:,dummy) = dev_1to50_cyl(:,dummy)/base_cyl(:,dummy);
MRE_1to100_cyl(:,dummy)= dev_1to100_cyl(:,dummy)/base_cyl(:,dummy);
MRE_1to250_cyl(:,dummy)= dev_1to250_cyl(:,dummy)/base_cyl(:,dummy);
MRE_1to500_cyl(:,dummy)= dev_1to500_cyl(:,dummy)/base_cyl(:,dummy);
    %%% Wake
MRE_1to1_wake(:,dummy)  = dev_1to1_wake(:,dummy)/base_wake(:,dummy);
MRE_1to5_wake(:,dummy)  = dev_1to5_wake(:,dummy)/base_wake(:,dummy);
MRE_1to10_wake(:,dummy) = dev_1to10_wake(:,dummy)/base_wake(:,dummy);
MRE_1to50_wake(:,dummy) = dev_1to50_wake(:,dummy)/base_wake(:,dummy);
MRE_1to100_wake(:,dummy)= dev_1to100_wake(:,dummy)/base_wake(:,dummy);
MRE_1to250_wake(:,dummy)= dev_1to250_wake(:,dummy)/base_wake(:,dummy);
MRE_1to500_wake(:,dummy)= dev_1to500_wake(:,dummy)/base_wake(:,dummy);
end

% temp_mean_fullfow=mean(temp_MRE_fullflow*100);
%Mean value of MRE Wake
Mean_MRE1to1_wake=mean(MRE_1to1_wake*100);
Mean_MRE1to10_wake=mean(MRE_1to10_wake*100);
Mean_MRE1to100_wake=mean(MRE_1to100_wake*100);
Mean_MRE1to5_wake=mean(MRE_1to5_wake*100);
Mean_MRE1to50_wake=mean(MRE_1to50_wake*100);
Mean_MRE1to250_wake=mean(MRE_1to250_wake*100);
Mean_MRE1to500_wake=mean(MRE_1to500_wake*100);
Mean_MRE_wake=[Mean_MRE1to1_wake Mean_MRE1to5_wake Mean_MRE1to10_wake ... 
    Mean_MRE1to50_wake Mean_MRE1to100_wake ... 
    Mean_MRE1to250_wake Mean_MRE1to500_wake];

%Mean value of MRE Inlet
Mean_MRE1to1_inlet=mean(MRE_1to1_inlet*100);
Mean_MRE1to10_inlet=mean(MRE_1to10_inlet*100);
Mean_MRE1to100_inlet=mean(MRE_1to100_inlet*100);
Mean_MRE1to5_inlet=mean(MRE_1to5_inlet*100);
Mean_MRE1to50_inlet=mean(MRE_1to50_inlet*100);
Mean_MRE1to250_inlet=mean(MRE_1to250_inlet*100);
Mean_MRE1to500_inlet=mean(MRE_1to500_inlet*100);
Mean_MRE_inlet=[Mean_MRE1to1_inlet Mean_MRE1to5_inlet Mean_MRE1to10_inlet ... 
    Mean_MRE1to50_inlet Mean_MRE1to100_inlet ... 
    Mean_MRE1to250_inlet Mean_MRE1to500_inlet];
%Mean value of MRE cylinder
Mean_MRE1to1_cyl=mean(MRE_1to1_cyl*100);
Mean_MRE1to10_cyl=mean(MRE_1to10_cyl*100);
Mean_MRE1to100_cyl=mean(MRE_1to100_cyl*100);
Mean_MRE1to5_cyl=mean(MRE_1to5_cyl*100);
Mean_MRE1to50_cyl=mean(MRE_1to50_cyl*100);
Mean_MRE1to250_cyl=mean(MRE_1to250_cyl*100);
Mean_MRE1to500_cyl=mean(MRE_1to500_cyl*100);
Mean_MRE_cyl=[Mean_MRE1to1_cyl Mean_MRE1to5_cyl Mean_MRE1to10_cyl ... 
    Mean_MRE1to50_cyl Mean_MRE1to100_cyl ... 
    Mean_MRE1to250_cyl Mean_MRE1to500_cyl];
%% Figures MRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lwdth=1.25;
saveFig = true;
Extension =false;   %true = PDF , false = PNG
Loop=false;
%%% Temporal mean of MRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScriptAspectFigures_Pre;
xaxis=[1 5 10 50 100 250 500];
hold on
plot(xaxis,Mean_MRE_inlet,'-o','LineWidth',Lwdth)
plot(xaxis,Mean_MRE_cyl,'-o','LineWidth',Lwdth)
plot(xaxis,Mean_MRE_wake,'-o','LineWidth',Lwdth)
xlabel('Total amount of modes'),ylabel('Temporal mean of MRE [%]')
legend('Domain inlet area','Domain cylinder area', 'Domain outlet area','Location','northeast')
temptitle=('Temporal mean of MRE for different domains');
hold off
figName = 'MRE';
ScriptAspectFigures_Post;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% wake area MRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% for some reason this one is giving errors. you can do this one manually
%%% wihtout problems
ScriptAspectFigures_Pre;
hold on
%  plot(Time(1,:),MRE_1to1_wake(1,:)*100,'-','LineWidth',Lwdth)
% plot(Time(1,:),MRE_1to5_wake(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to10_wake(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to50_wake(1,:)*100,'-','LineWidth',Lwdth)
% plot(Time(1,:),MRE_1to100_wake(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to250_wake(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to500_wake(1,:)*100,'-','LineWidth',Lwdth)
xlabel('Time (s)'),  ylabel('MRE (%)')
legend('Mode 1 to 10','Mode 1 to 50','Mode 1 to 250','Mode 1 to 500','Mode 1 to 100' ...
    ,'Mode 1 to 250','Mode 1 to 500','Location', 'northeast');
temptitle=('MRE for outlet domain');
xlim([Time(1,1) Time(1,end)]), ylim([0 110])
figName = 'MRE_Outlet';
ScriptAspectFigures_Post;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Upperduct /CYL area MRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScriptAspectFigures_Pre;
hold on
% plot(Time(1,:),MRE_1to1_cyl(1,:)*100,'-','LineWidth',Lwdth)
% plot(Time(1,:),MRE_1to5_cyl(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to10_cyl(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to50_cyl(1,:)*100,'-','LineWidth',Lwdth)
% plot(Time(1,:),MRE_1to100_cyl(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to250_cyl(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to500_cyl(1,:)*100,'-','LineWidth',Lwdth)
xlabel('Time (s)'),  ylabel('MRE (%)')
legend('Mode 1 to 10','Mode 1 to 50','Mode 1 to 250','Mode 1 to 500','Mode 1 to 100' ...
    ,'Mode 1 to 250','Mode 1 to 500','Location', 'northeast');
temptitle=('MRE for cylinder domain');
xlim([Time(1,1) Time(1,end)]), ylim([0 110])
figName = 'MRE_Cylinder';
ScriptAspectFigures_Post;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Lower duct /inlet area MRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ScriptAspectFigures_Pre;
hold on
% plot(Time(1,:),MRE_1to1_inlet(1,:)*100,'-','LineWidth',Lwdth)
% plot(Time(1,:),MRE_1to5_inlet(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to10_inlet(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to50_inlet(1,:)*100,'-','LineWidth',Lwdth)
% plot(Time(1,:),MRE_1to100_inlet(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to250_inlet(1,:)*100,'-','LineWidth',Lwdth)
plot(Time(1,:),MRE_1to500_inlet(1,:)*100,'-','LineWidth',Lwdth)
xlabel('Time (s)'),  ylabel('MRE (%)')
legend('Mode 1 to 10','Mode 1 to 50','Mode 1 to 250','Mode 1 to 500','Mode 1 to 100' ...
    ,'Mode 1 to 250','Mode 1 to 500','Location', 'northeast');
temptitle=('MRE for inlet domain');
xlim([Time(1,1) Time(1,end)]), ylim([0 110])
figName = 'MRE_Inlet';
ScriptAspectFigures_Post;
close all
%% Temporal coefficients visualization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lwdth=1.25;
% saveFig = true;
% Extension =true;   %true = PDF , false = PNG
% Loop=true; 
% 
% for temploop=1:50
% ScriptAspectFigures_Pre;
% hold on
% plot(Time(1,:),Coefficients(:,temploop),'-','LineWidth',Lwdth);
% % plot(Coefficients_ip(:,temploop),'-','LineWidth',Lwdth);
% % legend('Case 3','Interpolated Case','Location','northeast')
% xlabel('Time (s)'),ylabel('Magnitude')
% xlim([Time(1,1)-dTime Time(1,end)+dTime])
% ScriptAspectFigures_Post;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% MODES - comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FtSize            = 8; %fontsize
% % 
% for temploop=1:50
%  
% trisurf(tri,X,Y,Modes(:,temploop)')
% caxis([-0.05 0.05]), colormap(jet), 
% 
% hcb=colorbar;
% shading interp, view(2), axis equal tight
% 
% % ylim([-0.1 0.1])
% xlabel('X (m)'),ylabel('Y (m)')
% ylabel(hcb,'Magnitude')
% title(sprintf('Mode %d, Case 2',temploop))
% 
% set(gca,'FontSize',FtSize);
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 15 4]);
% % set(gcf,'Renderer','Painters');        % vector format
% set(gcf,'RendererMode','manual');      % so that the Renderer is manual
% viscircles([0.04 0],0.00125,'Color','k')
% % plot(0.04, 0, '.k', 'MarkerSize',100)
% print('-dpng','-r400',sprintf('image%d.png',temploop));
% drawnow
% end
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CENTERLINE PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X_cl = X(Centerline);
% Flow_field_avg=mean(Flow_field);
% Flow_field_avg_lim1to5=mean(Flow_field_lim1);
% Flow_field_avg_lim1to25=mean(Flow_field_lim1to3);

% % time plots
% a=[];
% for t1=1:100
% Pressure_cl = Pressure(Centerline,t1) - mean(Pressure(Centerline,t1));
% E=round((sum(Pressure_cl_pod(:,t1) ) - sum(Pressure_cl) )/ sum(Pressure_cl) *100);
% 
% plot(X_cl,Pressure_cl,'bo',X_cl,Pressure_cl_pod(:,t1),'r*');
% grid on, grid minor,
% xlabel('X-coordinate'), ylabel('Pressure [Pa]')
% ylim([-1500 1500])
% 
% textLabel = sprintf('Relative error= %d%', E );
% dim = [0.15 0.2 0 0];
% delete(a)
% a=annotation('textbox',dim,'String',textLabel,'FitBoxToText','on');
% legend('LES','POD')
% 
% drawnow
% pause(0.2)
% end

% averaged plots
% E2= sum(Flow_field_avg) - sum(P_avg) / sum(P_avg)         *100 ; %% Error
% E3= sum(Flow_field_avg_lim1to5) - sum(P_avg) / sum(P_avg) *100 ; %% Error
% E4= sum(Flow_field_avg_lim1to25) - sum(P_avg) / sum(P_avg)*100 ; %% Error

% figure
% plot(X_cl,P_avg(Centerline),'bo',X_cl,Flow_field_avg(Centerline),'rx')
% grid on, grid minor
% xlabel('X-coor [m]'),ylabel('Pressure [Pa]')
% legend('LES data','POD, mode 0 and 1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pressure plot StarCCM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
FlucPressure=Pressure;
% -mean(Pressure,2);
for p = 1:size(Pressure,2)
Time_sim = Time(1,1) + p*dTime;    
Z= FlucPressure(:,p);
trisurf(tri,X,Y,Z');
viscircles([0.04 0],0.00125,'Color','k')
caxis([-10200 -9500]);

% caxis([-200 200]) 
colormap(jet), hcb=colorbar;
shading interp, view(2), axis equal tight

xlim([0 0.1])
xlabel('X (m)'),ylabel('Y (m)')
ylabel(hcb,'Pressure (Pa)')
title(sprintf('Pressure, (Time = %d s)',Time_sim))

set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 18 5.5]);
print('-dpng','-r0',sprintf('image%d.png',p));
end


% 
% Pressure_point_1801_avg=mean(Pressure_point_1801);
% figure
% plot(Time(1,:),Pressure_point_1801(1,:),'-o')
% hline=refline([0 Pressure_point_1801_avg]); hline.Color = 'r';
% grid on; grid minor
% ylim([-10100 -9500]),xlim([0.1 0.1271])
% xlabel('Time [s]'), ylabel('Pressure [Pa]')
% 
% Pressure_point_1801_POD_avg=mean(Pressure_point_1801_POD);
% figure
% plot(Time(1,:),Pressure_point_1801_POD(1,:),'-or')
% hline=refline([0 Pressure_point_1801_POD_avg]); hline.Color = 'g';
% grid on; grid minor
% xlim([0.1 0.1271])
% xlabel('Time [s]'), ylabel('Pressure [Pa]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Temp Probe point outside centerline -- REFLECTION research
% define new p.p
tempidx=find(abs(X-0.05)<1e-3 & abs(Y-0.1)<1e-3);
PP_od=Pressure(tempidx,:);

tempidx2=find(abs(X-0.65)<3e-3 & abs(Y-0.1)<1e-2); 
PP_od2=Pressure(tempidx2,:);

% tempidx=find(abs(X+0.25)<3e-3 & abs(Y-0.1)<1e-2 );
% PP_od=Pressure(tempidx(1),:);
% 
% tempidx2=find(abs(X+0.25)<3e-3 & abs(Y+0.1)<1e-2); 
% PP_od2=Pressure(tempidx2(1),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PP_od_fft=fft(PP_od);
% PP_od_fft2=fft(PP_od2);
% 
% Pds_PPod = abs(PP_od_fft/Lts);                     % double-sided amplitude spectrum
% Pss_PPod = Pds_PPod(1:Lts/2+1);                   % single-sided amplitude spectrum
% Pss_PPod(2:end-1) = 2*Pss_PPod(2:end-1);          % frequency shift
% 
% Pds_PPod2 = abs(PP_od_fft2/Lts);                  % double-sided amplitude spectrum
% Pss_PPod2 = Pds_PPod2(1:Lts/2+1);                 % single-sided amplitude spectrum
% Pss_PPod2(2:end-1) = 2*Pss_PPod2(2:end-1); 

% Extension=false;
% figure
% ScriptAspectFigures_Pre;
% hold on
% plot(f,Pss_PPod,'-','LineWidth',Lwdth)
% plot(f,Pss_PPod2,'-','LineWidth',Lwdth)
% xlabel('Frequency (Hz)'),ylabel('Pressure (Pa)')...
%     ,xlim([0 100]),ylim([0 100])
% legend('P.P upper duct','P.P lower duct','Location','northeast')
% temptitle=('Comparison between probe points');
% hold off
% figName = 'Dummytemp';
% ScriptAspectFigures_Post;
% 
% figure
% hold on
% triplot(tri,X,Y)
% % viscircles([X(PP050(1)) Y(PP050(1))],0.005,'Color','r')
% viscircles([X(tempidx(1)) Y(tempidx(1))],0.005,'Color','r')
% % viscircles([X(PP650(1)) Y(PP650(1))],0.005,'Color','g')
% viscircles([X(tempidx2(1)) Y(tempidx2(1))],0.005,'Color','g')
% axis equal tight
% xlim([-0.29 0.05])
% xlabel('X (m)'); ylabel('Y (m)')
% title('Probe point comparison')
% % 
% Lwdth=1.25;
% saveFig = true;
% Extension =true;   %true = PDF , false = PNG
% Loop=false;
% 
% ScriptAspectFigures_Pre;
% hold on
% plot(Time(1,:),Pressure_PP050,'-','LineWidth',Lwdth)
% plot(Time(1,:),PP_od,'-','LineWidth',Lwdth)
% % plot(Time(1,:),PP_od2,'-','LineWidth',Lwdth)
% xlim([Time(1,1)-dTime Time(1,end)+dTime]);ylim([-200 200])
% xlabel('Time (s)'),ylabel('Pressure (Pa)')
% legend('P.P @ X=-250 mm and Y=+100 mm','P.P @ X=-250 mm and Y=-100 mm','Location','northeast')
% temptitle=('Pressure comparison');
% figName = 'PressureComparisonWake_Domain_X50mm';
% ScriptAspectFigures_Post;
% 
% ScriptAspectFigures_Pre;
% hold on
% plot(Time(1,:),Pressure_PP650,'-','LineWidth',Lwdth)
% plot(Time(1,:),PP_od2,'-','LineWidth',Lwdth)
% xlim([Time(1,1)-dTime Time(1,end)+dTime]);ylim([-200 200])
% xlabel('Time (s)'),ylabel('Pressure (Pa)')
% legend('P.P @ X=650 mm and Y=0 mm','P.P @ X=650 mm and Y=100 mm','Location','northeast')
% temptitle=('Pressure comparison');
% figName = 'PressureComparisonWake_Domain_X650mm';
ScriptAspectFigures_Post;
%% Temp Interpolated flow field vs LES flow field - MRE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error plots - comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dTime=Time(1,2)-Time(1,1);
tallocation=size(Time,2);
ngrid=size(X_coor,1);
%%%% define wake area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=X_coor(:,1);
Y=Y_coor(:,1);
LLimit=0.035; ULimit=0.05;
%%% find indexes for mixing layer problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yud_upper_ul=0.14; Yud_upper_ll=0.06;                   % upper duct limits
idx_upperduct=find(Y>Yud_upper_ll & Y<Yud_upper_ul);
X_upper=X(idx_upperduct);Y_upper=Y(idx_upperduct);

Yld_lower_ll=-0.14; Yld_upper_ll=-0.06;                 % lower duct limits
idx_lowerduct=find(Y>Yld_lower_ll & Y<Yld_upper_ll);
X_lower=X(idx_lowerduct);Y_lower=Y(idx_lowerduct);

idx_wake=find(Y<Yud_upper_ll & Y>Yld_upper_ll & X>-0.3);
X_wake=X(idx_wake);Y_wake=Y(idx_wake);

idx_wake=idx_wake;idx_cyl=idx_upperduct;idx_inlet=idx_lowerduct;
%%%% We are interested in the MRE of the fluctuating field %%%%%%%%%%%%%%%%
Pressure_fluc=Pressure-mean(Pressure,2);

temp_ip_wake=abs(Pressure_fluc(idx_wake,:) - Flow_field_ip(:,idx_wake)');
dev_ip_wake=sum(temp_ip_wake,1);

temp_ip_cyl=abs(Pressure_fluc(idx_cyl,:) - Flow_field_ip(:,idx_cyl)');
dev_ip_cyl=sum(temp_ip_cyl,1);

temp_ip_inlet=abs(Pressure_fluc(idx_inlet,:) - Flow_field_ip(:,idx_inlet)');
dev_ip_inlet=sum(temp_ip_inlet,1);

%%% Base for the three domains
base_inlet=sum(abs(Pressure_fluc(idx_inlet,:)));
base_cyl=  sum(abs(Pressure_fluc(idx_cyl,:)));
base_wake= sum(abs(Pressure_fluc(idx_wake,:)));

for dummy= 1:tallocation
MRE_ip_inlet(:,dummy)=dev_ip_inlet(:,dummy)/base_inlet(:,dummy);
MRE_ip_cyl(:,dummy)  =dev_ip_cyl(:,dummy)/base_cyl(:,dummy);
MRE_ip_wake(:,dummy) =dev_ip_wake(:,dummy)/base_wake(:,dummy);
end

MeanMRE_ip_inlet=mean(MRE_ip_inlet*100);
MeanMRE_ip_cyl=mean(MRE_ip_cyl*100);
MeanMRE_ip_wake=mean(MRE_ip_wake*100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lwdth=1.25;
saveFig = true;
Extension =false;   %true = PDF , false = PNG
Loop=false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% wake area MRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ScriptAspectFigures_Pre;
hold on
plot(MRE_ip_wake(1,:)*100,'-','LineWidth',Lwdth)
xlabel('Time-step'),  ylabel('MRE (%)')
% legend('Mode 1 to 1','Mode 1 to 5','Mode 1 to 10','Mode 1 to 50','Mode 1 to 100' ...
%     ,'Mode 1 to 250','Mode 1 to 500','Location', 'northeast');
temptitle=('MRE for wake domain');
% xlim([Time(1,1) Time(1,end)]), ylim([0 110])
figName = 'MRE_wake';
ScriptAspectFigures_Post;
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Upperduct /CYL area MRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScriptAspectFigures_Pre;
hold on
plot(MRE_ip_cyl(1,:)*100,'-','LineWidth',Lwdth)
xlabel('Time-step'),  ylabel('MRE (%)')
% legend('Mode 1 to 1','Mode 1 to 5','Mode 1 to 10','Mode 1 to 50','Mode 1 to 100' ...
%     ,'Mode 1 to 250','Mode 1 to 500','Location', 'northeast');
temptitle=('MRE for upper duct domain');
% xlim([Time(1,1) Time(1,end)]), ylim([0 110])
figName = 'MRE_upperduct';
ScriptAspectFigures_Post;
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Lower duct /inlet area MRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ScriptAspectFigures_Pre;
hold on
plot(MRE_ip_inlet(1,:)*100,'-','LineWidth',Lwdth)
xlabel('Time-step'),  ylabel('MRE (%)')
% legend('Mode 1 to 1','Mode 1 to 5','Mode 1 to 10','Mode 1 to 50','Mode 1 to 100' ...
%     ,'Mode 1 to 250','Mode 1 to 500','Location', 'northeast');
temptitle=('MRE for lower duct domain');
% xlim([Time(1,1) Time(1,end)]), ylim([0 110])
figName = 'MRE_lowerduct';
ScriptAspectFigures_Post;
close all
%% Figures Report 
% cd('F:\MThesis\UNIX\Desktop\CFD\2D_KarmanDOE\DOE\FiguresReport\Param1_50/Vorticity')
cd('F:\MThesis\UNIX\Desktop\CFD\2D_KarmanDOE\DOE\FiguresReport\SpatialModesCompare\Param2_00')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm=1000;
% ScriptAspectFigures_Pre_Report;
% trisurf(tri, X*mm,Y*mm, mu') 
% view(2); axis equal tight; hcb=colorbar; shading interp; colormap jet
% xlim([0.03*mm 0.1*mm]); ylim([-0.01*mm 0.01*mm]);
% xlabel('X (mm)'); ylabel('Y (mm)');
% caxis([-1.1e4 -0.9e4])
% % caxis([0 8e4])
% viscircles([0.04*mm 0],0.00125*mm,'Color','k');
% title('Mean flow-field')
% ylabel(hcb, 'Pressure (Pa)')
% hcb.Label.FontSize=10;
% 
% figName = 'MeanField';
% Loop=false;
% saveFig =true;    %true = safe, false = don't safe
% Extension = false;   %true = PDF , false = PNG
% ScriptAspectFigures_Post_Report;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=13;
Modes(:,:)=Param_POD{1,i}{1,1};
for temploop=10:25
ScriptAspectFigures_Pre_Report;
trisurf(tri, X*mm,Y*mm, Modes(:,temploop)) 
view(2); axis equal tight; hcb=colorbar; shading interp; colormap jet
xlim([0.03*mm 0.1*mm]); ylim([-0.01*mm 0.01*mm]);
caxis([-0.02 0.02])
xlabel('X (mm)'); ylabel('Y (mm)');
viscircles([0.04*mm 0],0.00125*mm,'Color','k');



Loop=true;

ScriptAspectFigures_Post_Report;
end

%% Temporal plots figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('F:\MThesis\UNIX\Desktop\CFD\2D_KarmanDOE\DOE\FiguresReport\TemporalPlots2\Interpolated\Param=1_50_Ip')

Lwdth=1.25;
saveFig = true;
Extension =true;   %true = PDF , false = PNG
Loop=false;
dTime=Time(1,2)-Time(1,1);

for temploop=1:25
ScriptAspectFigures_Pre_TemporalFigures;
hold on
plot(Time(1,:),Coefficients_c1(:,temploop),'-','LineWidth',Lwdth);
% plot(Time(1,:),Coefficients_c1_5(:,temploop),'-','LineWidth',Lwdth);
plot(Time(1,:),Coefficients_ip(:,temploop),'-','Color',[0.8500 0.3250 0.0980],'LineWidth',Lwdth);
plot(Time(1,:),Coefficients_c2(:,temploop),'-','Color',[ 0.3010    0.7450    0.9330],'LineWidth',Lwdth);
xlim([0.15 0.16])
% plot(Coefficients_ip(:,temploop),'-','LineWidth',Lwdth);
legend('Training point 1','Design point 8','Training point 13','Location','northeast')
xlabel('Time (s)'),ylabel('Magnitude')
% xlim([Time(1,1)-dTime Time(1,end)+dTime])
figName=(sprintf('Temporal_mode%d', temploop));
ScriptAspectFigures_Post_TemporalFigures;
end




