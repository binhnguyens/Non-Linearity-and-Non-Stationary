%%%%%%%%% Non-linearity test %%%%%%%%%
% Non-linear if p < 0.05
% Linear if p > 0.05
clear; clc;
n = 9;
path = '../Data/';


switch n 

     case 1
         
        %PPG Signal
        
        load (strcat (path, 'PPG_signal.mat'));
       
        signal = normalization (PPG (1:15000));
        fs = 2000;
        
     case 2
        
         % AR Synthetic Signal
        synthetic_AR = zeros (1,1000);
        synthetic_AR(1) = 1000;
        for i = 2:1000
            synthetic_AR (i) = 0.5*synthetic_AR(i-1) ;
        end

        signal = transpose (normalization (synthetic_AR (1:100)));
        fs = 1000;
        
    case 3
        % Linear Signal
        linear_line = 1:1000;
        
        signal = transpose (normalization (linear_line));
        fs = 1;
        
    case 4
        % MN Signal
        filename = 'test_ES2004c_1980_1990_xxxxx_4_0.wav';
        y = audioread (strcat (path , filename));
        
        
        signal = normalization (y(1:20000));
        fs = 16000;
        
    case 5
        % AR Signal
        filename = 'AVPEPUDEAC0045a1.wav';
        y = audioread (strcat (path , filename));
        
        signal = normalization (y(1:20000));
        fs = 44000;
    
    case 6
        % Sinusoid signal
           fs = 1000; 
           
           dt = 1/fs; 
           StopTime = 0.25;  
           t = (0:dt:StopTime-dt)';   
           
           Fc = 60;                   
           x = cos(2*pi*Fc*t);
           
           signal = normalization (x);
           
     case 7
         
        %ECG Signal
        
        load (strcat (path, 'ECG.mat'));
       
        signal = normalization (ECG (1:15000));
        fs = 2000;
        
        
     case 8
         
        %RI Signal
        
        load (strcat (path, 'RI.mat'));
       
        signal = normalization (RI (1:15000));
        fs = 2000;
                
     case 9
         
        %GSR Signal
        
        load (strcat (path, 'GSR.mat'));
       
        signal = normalization (GSR (1:15000));
        fs = 2000;
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

input = 1:length(signal);
z = iddata(signal,transpose (input), 1/fs);
order = [0 1 0];
isnlarx (z,order)
