  
clear; close all; clc;

for par_num = 1:1

    %% Add path file to CSV or TXT files
    % Main file
    path = ('/Users/binhnguyen/Desktop/Desktop/Digital Mental Health/2. Data and Analysis/Biopac Data/');

    % Participant files and Event markers


    switch par_num
        case 1
            filename1 = ('Session 1 - May 31/2021-05-31 Participant 1.csv');
            filename2 = ('Session 1 - May 31/event markers p1.csv');

        case 2
            filename1 = ('Session 1 - May 31/2021-05-31 Participant 2.csv');
            filename2 = ('Session 1 - May 31/event markers p2.csv');

        case 3
            filename1 = ('Session 1 - May 31/2021-05-31 Participant 3.csv');
            filename2 = ('Session 1 - May 31/event markers p3.csv');

        case 4
            filename1 = ('Session 2 and 3 - June 8 and 10/2021-06-08 - Participant 04.csv');
            filename2 = ('Session 2 and 3 - June 8 and 10/event markers p4.csv');

        case 5
            filename1 = ('Session 2 and 3 - June 8 and 10/2021-06-08 - Participant 05.csv');
            filename2 = ('Session 2 and 3 - June 8 and 10/event markers p5.csv');

        case 6
            filename1 = ('Session 2 and 3 - June 8 and 10/2021-06-10 - Participant 06.csv');
            filename2 = ('Session 2 and 3 - June 8 and 10/event markers p6.csv');

        case 7
            filename1 = ('Session 4 - June 15/2021-06-15 - Participant 07.csv');
            filename2 = ('Session 4 - June 15/event markers p7.csv');    

        case 8
            filename1 = ('Session 5 - June 25/2021-06-25 - Participant 08.csv');
            filename2 = ('Session 5 - June 25/event markers p8.csv');

        case 9
            filename1 = ('Session 6-8/session6-8 - Participant 09.csv');
            filename2 = ('Session 6-8/event markers p9.csv');

        case 10
            filename1 = ('Session 6-8/session6-8 - Participant 10.csv');
            filename2 = ('Session 6-8/event markers p10.csv');

        case 11
            filename1 = ('Session 6-8/session6-8 - Participant 11.csv');
            filename2 = ('Session 6-8/event markers p11.csv');

        case 12
            filename1 = ('Session 6-8/session6-8 - Participant 12.csv');
            filename2 = ('Session 6-8/event markers p12.csv');

        case 13
            filename1 = ('Session 6-8/session6-8 - Participant 13.csv');
            filename2 = ('Session 6-8/event markers p13.csv');

        case 14
            filename1 = ('Session 6-8/session6-8 - Participant 14.csv');
            filename2 = ('Session 6-8/event markers p14.csv');

        case 15
            filename1 = ('Session 6-8/session6-8 - Participant 15.csv');
            filename2 = ('Session 6-8/event markers p15.csv');

        otherwise
            fprintf ("False response");

    end

    file = strcat ( path , filename1 );

    event_marker_file = strcat ( path , filename2 );


    %% Open files
    df = readtable (file);

    % Retrieve the event markers
    df_event_markers = readtable (event_marker_file);
    df_EM = event_marker_updated (df_event_markers);

    % Remove first row
    df(1,:) = [];

    %% Pre-processing

    % Participant 7 - data removal at the end 

    if (par_num == 7)

        df(6760420:end,:) = [];

    end
    
    % Normalization
    
    n_c = size (df,2);

    for i = 1:n_c
        df{:,i} = normalization (df{:,i});
    end

    %% Assign value
    RI = normalization (df.CH1);
    PPG = normalization (df.CH2);
    ECG = normalization (df.CH13);
    GSR = normalization (df.CH14);
    
end

%% Non-linear party

% Steps
% 1. Separate the segments into 5 because it has to be normalized
% 2. Since the first order is very stationary, we can analyze the entire
% signal
% 3. The signal is only stationary in the first order
% 4. We can only analyze signals or segments that are stationary 
% 5. We need to reduce the segments to make it segment (make the 10^-3 as a
% STD of means)
% 6. Finding the moment that the signal goes from stationary to
% non-stationary then that can be considered as NON-STATIONARY (Region of
% interest)
% 7. We can also look into high-order nonstationary (e.g first order, second
% order, and higher such as kurtosis or skewness)

X = ECG (1:5500000); 
X_reshaped = reshape(X,[],5);

% Mean value of segments
% Signal is stationary bc mean is consistent 
X_mean_of_segments = mean (X_reshaped);
X_STD_of_Means = std (X_mean_of_segments);

% Variance value of segments
% Signal is stable/stationary
X_var_of_segments = var (X_reshaped);
X_STD_of_VAR = std (X_var_of_segments);


fprintf ('First order analysis\n');
disp (X_mean_of_segments);
fprintf ('X STD of means is %d \n', X_STD_of_Means);

disp ('-------------------------------------------');

fprintf ('Second order analysis\n');
disp (X_var_of_segments);
fprintf ('X STD of means is %d \n', X_STD_of_VAR);


disp ('-------------------------------------------');
fprintf ('Duration of signal %.2d mins\n', length(X)/2000/60);


% BIG TAKEAWAY
% Since signals are stationary, then we do not need to use complex models
% since it is stationary and we can analyze the entire signal

