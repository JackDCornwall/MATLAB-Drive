%%
% Author: Robin Smith    Date: 12/11/2015 %%%%

% Polynomial fitting program %%%%
%
% A code that reads in a xy column file   %%%%
% and fits a polynomial to the data.      %%%%
% Fit is weighted by y error bars         %%%%

% READ IN DATA %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

filename=uigetfile;

Dat=dlmread(filename);

X = Dat(:,1);

Y = Dat(:,2);

Err = Dat(:,3);

Size=size(X);

NumPoints=Size(1);

XlabIn=input('Enter X axis Label (use inverted commas): ');
Xlab=XlabIn;

YlabIn=input('Enter Y axis Label (use inverted commas): ');
Ylab=YlabIn;

% FIT DATA %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Weight=1./Err;

Wfunc=@(p,input) Weight.*LinearPoly(p,input);

WY=Weight.*Y;

init=zeros(1,2);

% Approximate initial fit parameters
init(1)=(Y(NumPoints)-Y(1))/ (X(NumPoints)-X(1));

init(2)=Y(1)-init(1)*X(1);

options=optimset('lsqcurvefit');
options.MaxFunEvals=50000;
options.MaxIter=5000;
options.TolFun=1e-9;

[NewParameters,error,residual,exitflag,output,lambda,Jac]=lsqcurvefit(Wfunc,init,X,WY);

newY=LinearPoly(NewParameters,X);

% PLOT DATA AND FIT LINE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(X,newY);
xlabel(Xlab,'FontSize',14);
ylabel(Ylab,'FontSize',14);
%title('Resolution vs 9Be Excitation: Monte Carlo Simulations','FontSize',14);
hold;

errorbar(X,Y,Err,'.k');

% CALCULATE ERRORS ON THE FIT PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resid=newY-Y;
[Q,R] = qr(Jac,0);
mse = sum(abs(resid).^2)/(size(Jac,1)-size(Jac,2));
Rinv = inv(R);
Sigma = Rinv*Rinv'*mse;
se = sqrt(diag(Sigma));

Gradient=['Gradient ' num2str(NewParameters(1)) ' +/- ' num2str(se(1))]
Intercept=['Intercept ' num2str(NewParameters(2)) ' +/- ' num2str(se(2))]

% CALCULATE THE GOODNESS OF FIT CHI2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chi2 = 0;

for i=1:length(X)
    chi2 = chi2 + (abs(resid(i))./abs(Err(i).^2));
end

ChiSquared=['Chi2 = ' num2str(chi2)]
ChiSquared=['Chi2 per degree of freedom = ' num2str(chi2/(length(X) - 2))]