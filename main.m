%%%%%%%%% Non-linearity test %%%%%%%%%
% Non-linear if p < 0.05
% Linear if p > 0.05
clear; clc;
n = 5;

switch n 

     case 1
         
        %PPG Signal
        
        load ('PPG_signal.mat');
        signal = normalization (PPG (1:15000));
        
        
     case 2
        
         % AR Signal
        synthetic_AR = zeros (1,1000);
        synthetic_AR(1) = 1000;
        for i = 2:1000
            synthetic_AR (i) = 0.5*synthetic_AR(i-1) ;
        end

        signal = normalization (synthetic_AR (1:100));
        
    case 3
        % Linear Signal
        linear_line = 1:100;
        signal = normalization (linear_line);
        
    case 4
        % MN Signal
        path = '/Users/binhnguyen/Downloads/';
        filename = 'test_ES2004c_1980_1990_xxxxx_4_0.wav';
        y = audioread (strcat (path , filename));
        
        
        signal = normalization (y(1:20000));
        
    case 5
        %AR Signal
        path = '/Users/binhnguyen/Downloads/';
        filename = 'AVPEPUDEAC0045a1.wav';
        y = audioread (strcat (path , filename));
        
        signal = normalization (y(1:20000));

end

%% Fast BDS test 

% Set values
M = 7;
EPS = 0.5;
RAMSIZE = 300;

% BDS test for Indepdence (w statistical value)
[W, SIG, C, C1, K] = bds (signal, 7, EPS, 0, RAMSIZE);

% BDS for small test
SIG = bdssig (W, length (signal), 7, EPS);

fprintf ('This is the p-value test for BDS\n');
disp (SIG);


%% isnarlx test
% Larger values (>2) indicate that a significant nonlinearity was detected.
% 
% Smaller values (<0.5) indicate that any error unexplained by the linear model is mostly noise. That is, no significant nonlinearity was detected.
% 
% Values close to 1 indicate that the nonlinearity detection test is not reliable and that a weak nonlinearity may be present.



