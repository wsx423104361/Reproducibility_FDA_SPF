clc;
clear;

%% Load Data
Density_h1     = readmatrix('.\Data\Density_h1.csv');
LQD_h1         = readmatrix('.\Data\LQD_h1.csv');
UniDate        = datetime(1981,09,30):calquarters(1):datetime(2022,06,30);
UniDate.Format = 'yyyy_QQQ';
nDate          = length(UniDate);

Range1         = -5:(15/511):10;
Range2         = 0:(1/255):1;

%% Plot
figure
colormap('parula');
Sparse = 171:2:512; % Plot from 0 to 10 in the x-axis with sparse datapoints for better visualization of density curves
surf(UniDate, Range1(Sparse), Density_h1(Sparse,:))
xlabel('Date')
zlabel('Density')

figure
colormap('hot');
surf(UniDate, Range2, LQD_h1)
xlabel('Date')
zlabel('LQD')