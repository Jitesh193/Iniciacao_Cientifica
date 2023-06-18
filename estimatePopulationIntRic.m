function [theta,Q]=estimatePopulationIntRic(x,y,id,h,theta0)
% Solve the nonlinear mixed-effect problem
% using Stochastic Approximation Expectation-Maximization (SAEM)
[BETA,Q] = nlmefitsa(x,y,id,[],h,theta0,...
    'ErrorModel','constant',... % measurement additive noise only
    'ParamTransform',[1 1 1 1],... % parameters are log-normal 
    'ComputeStdErrors',false);%,...
%     'LogLikMethod','is');
theta=exp(BETA); %fixed-effect
Q=diag(Q); %covariance matrix of the random effect
end