function [fitresult, gof] = CAfit(XData, YData,Tpv)
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
ft = fittype('1+(4*(Tpv/(0.406*(1-S)^0.27))*(x-t))/((1+(x-t)^2)*(9+(x-t)^2))-0.08  ', 'independent', 'x', 'dependent', 'y','problem','Tpv' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.15 0];
%             zero Slope Treshold threshold2 slop2
opts.StartPoint = [0.25 0.5];
opts.Upper = [0.3 1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts,'problem',Tpv);




