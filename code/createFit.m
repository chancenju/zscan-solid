function [fitresult, gof] = createFit(XData, YData)
%CREATEFIT(XDATA,YDATA)
%  Create a fit.
%
%  Data for 'threshold fit' fit:
%      X Input : XData
%      Y Output: YData
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 20-Oct-2020 11:46:17 自动生成

% *heaviside(threshold3-x)+(slope3*x+slope3*threshold3)*heaviside(x-threshold3)
%% Fit: 'threshold fit'.
[xData, yData] = prepareCurveData( XData, YData );

% Set up fittype and options.
ft = fittype('1-q0/(2*sqrt(2))*(1/(1+x^2/z0^2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
%             zero Slope Treshold threshold2 slop2
opts.StartPoint = [0.8491 0.9340 ];
opts.Upper = [Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );




