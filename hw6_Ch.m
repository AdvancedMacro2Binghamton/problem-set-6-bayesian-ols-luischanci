%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% PhD in Economics                              %
% ECON634 Advanced Macroeconomics               %
% Fall 2017                                     %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data:
cd 'C:\Users\Chanci\Dropbox\PhD\III. Third Year\1. ECON-634 Advanced Macroeconomics\Hw\hw6'
data  = csvread('data1.csv');
[n k] = size(data);
Y     = data(:,1);
X     = [ones(n,1) data(:,2:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target D. and Other functions
L = @(T,V,E)(-(T/2)*log(2*pi)-(T/2)*log(V)-inv(2*V)*(E'*E));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. OLS
b  = (inv(X'*X)*(X'*Y));
s  = inv(n-k)*(Y-X*b)'*(Y-X*b);
Vb = s*diag(inv(X'*X));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Flat Prior
r      = 0.002; % I played with this number to have acceptance rate 20-25% 
Sigma  = r*[[bsxfun(@times,Vb,eye(k)) zeros(k,1)];[zeros(1,k) s]];

% 2.1. Firt, track accept-reject status
B      = [b;s];
acc0    = [0,0];
for i  = 1:1e4                          % MH routine; 
    [B,a] = MHstep_Ch(B,Sigma,Y,X,L,1); % Option 1: flat prior
    acc0   = acc0 + [a 1];                % track accept-reject status
end
    acc0(1)/acc0(2)*100      % Acceptance rate. IF it is low, increase r.

% 2.2. Second, MH routine (after the burn-in)
lag    = 1;
nsamp  = 1e6;              % number of samples to draw
Theta  = zeros(nsamp,k+1); % storage
acc1    = [0,0];
for i = 1:nsamp
    for j=1:lag
        [B,a] = MHstep_Ch(B,Sigma,Y,X,L,1); % Option 1: flat prior
        acc1 = acc1 + [a 1];
    end
    Theta(i,:) = B';
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.   Using the Prior for Beta - Edu

% 3.1. Firt, track accept-reject status
B      = [b;s];
acc2   = [0,0];
for i  = 1:1e4                          % MH routine; 
    [B,a] = MHstep_Ch(B,Sigma,Y,X,L,2); % Option 2: Prior
    acc2  = acc2 + [a 1];                % track accept-reject status
end
    acc2(1)/acc2(2)*100  % Acceptance rate. IF it is low, increase r.

% 3.2. Second, MH routine
Theta2 = zeros(nsamp,k+1);  % storage
acc3   = [0,0];
for i = 1:nsamp
    for j=1:lag
        [B,a] = MHstep_Ch(B,Sigma,Y,X,L,2); % Option 2: Prior
        acc3  = acc3 + [a 1];
    end
    Theta2(i,:) = B';
end
%%%%%%%%%%%%%%    Now we can Plot  %%%%%%%%%%%%%%
