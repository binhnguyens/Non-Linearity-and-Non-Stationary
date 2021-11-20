% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % Executable part of main function BDS.M ends here % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %
% The following sub-function is not actually used by the main function and only %
% included for the benefit of those who would like to implement the BDS test in a %
% language which is either incapable of or inefficient in handling bit-wise AND- %
% operations, or those who would like to cross-check the above computation. Deleting %
% the sub-function from the script will NOT result in any increase in performance. %
% %
% To use the function, save the remainder of this code in a file named BDSNOBIT.M. %
% %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function w = bdsnobit (series, maxdim, eps)
%BDSNOBIT BDS test for independence IMPLEMENTED WITHOUT USING BIT-WISE FUNCTIONS
%
% Only such comments which relate exclusively to this implementation of the test and
% which cannot be found in the main function are included below.
%
% Copyright (c) 14 April 1998 by Ludwig Kanzler
% Department of Economics, University of Oxford
% Postal: Christ Church, Oxford OX1 1DP, England
% E-mail: ludwig.kanzler@economics.oxford.ac.uk
% $ Revision: 1.3 $ $ Date: 30 April 1998 $
%%%%%%%%%%%%%%%%%%%%%% Check and transformation of input arguments %%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
eps = 1;
if nargin == 1
maxdim = 2;
elseif maxdim < 2
error('MAXDIM needs to be at least 2!');
end
end
epsilon = std(series)*eps;
series = series(:)'; % PAIRS is the total number of unique pairs which can be
n = length(series); % formed from all observations (note that while this is
pairs = sum(1:n-1); % just (N-1)*N/2, MATLAB computes SUM(1:N-1) twice as fast!)
%%%%%%%%%%%% Computation and storage of one-dimensional distance information %%%%%%%%%%%%
% Recall that in the implementation of the main function above, table C is stored in bit-
% representation. When this is not possible or desirable, the second best method is to use
% one continuous vector of unassigned 8-bit integers (called UINT8). This, however,
% requires version 5.1 or higher, and a similar option may not be available in other high-
% level languages. Implementation does not depend on the ability to use unassigned low-bit
% integers and would work equally with double-precision integers, but the memory
% requirements would, of course, be higher. Using UINT8's is still a rather inefficient
% way of storing zeros and ones, which in principle require only a single bit each. On the
% PC, MATLAB actually requires "only" around 5 bytes for each UNIT8.
b(1:pairs) = uint8(0);
for i = 1 : n-1
b(1+(i-1)*(n-1)-sum(0:i-2):i*(n-1)-sum(1:i-1)) = abs(series(i+1:n)-series(i))<=epsilon;
end
clear series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computation of parameter K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sums(1 : n) = 0;
for i = 1 : n
sums(i) = sum(b(i+(0 : i-2)*n - cumsum(1 : i-1)))... % sum over column I
+ 1 ... % diagonal element
+ sum(b(1+(i-1)*(n-1)-sum(1:i-2) : i*(n-1)-sum(1:i-1))); % sum over row I
end
k = (sum(sums.^2) + 2*n - 3*(2*sum(b)+n)) / n/(n-1)/(n-2);
%%%%%%%%%%%%%%%%%% Computation of one-dimensional correlation estimates %%%%%%%%%%%%%%%%%%
bitsum(1:maxdim) = sum(b(1+(maxdim-1)*(n-1)-sum(0:maxdim-2) : pairs));
for m = maxdim-1 : -1 : 1
bitsum(m) = bitsum(m+1) + sum(b(1+(m-1)*(n-1)-sum(0:m-2):m*(n-1)-sum(1:m-1)));
end
c1(maxdim:-1:1) = bitsum(maxdim:-1:1) ./ cumsum([sum(1:n-maxdim), n-maxdim+1 : n-1]);
%%%%%%%%%% Computation of correlation estimates and SIGMA for higher dimensions %%%%%%%%%%
for m = 2 : maxdim
% Indexing in vector space once again follows the rules set out above. Multiplication
% is done by moving up column by column into north-west direction, so counter I runs
% backwards in the below WHILE loop until the Mth column (from the left) is reached:
i = n;
while i - m
% Multiplication is not defined on UINT8 variables and translating the columns
% twice, once from UINT8 to DOUBLE integer and then back to UINT8, would be
% inefficient, so it is better to sum entries (this operation - undocumented by
% MATLAB - is defined, and even faster than the documented FIND function!) and
% compare them against the value 2:
b(i + (m-1 : i-2)*n - sum(1:m-1) - cumsum(m : i-1)) = ...
sum([ b(i + (m-1 : i-2)*n - sum(1:m-1) - cumsum(m : i-1)); ...
b(i-1 + (m-2 : i-3)*n - sum(1:m-2) - cumsum(m-1 : i-2)) ]) == 2;
% The sum over each column is computed immediately after that column has been
% updated. To store the column sums, the vector SUMS already used above for the row
% sums is recycled (this is more memory-efficient than clearing the above SUMS
% vector and defining a new vector of the column sums, because in the latter case,
% MATLAB's memory space will end up being fragmented by variables K and C added to
% the memory in the meantime!):
sums(i) = sum(b(i + (m-1 : i-2)*n - sum(1:m-1) - cumsum(m : i-1)));
i = i - 1;
end
c(m-1) = sum(sums(m+1:n)) / sum(1:n-m);
sigma(m-1) = 2*sqrt(k^m + 2*k.^(m-(1:m-1))*(c1(1).^(2*(1:m-1)))'... % could use above
+ (m-1)^2*c1(1)^(2*m) - m^2*k*c1(1)^(2*m-2)); % inter-power subend
% functions instead
%%%%%%%%%%%%%%% Computation of the BDS statistic and level of significance %%%%%%%%%%%%%%%
w = sqrt(n-(2:maxdim)+1) .* (c-c1(2:maxdim).^(2:maxdim)) ./ sigma; % or use sub-functions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % Sub-function BDSNOBIT.M ends here % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% REFERENCES:
%
% Brock, William, Davis Dechert & José Scheinkman (1987), "A Test for Independence
% Based on the Correlation Dimension", University of Wisconsin-Madison, Social
% Science Research Working Paper, no. 8762
%
% Brock, William, Davis Dechert, José Scheinkman & Blake LeBaron (1996), "A test for
% independence based on the correlation dimension", Econometric Reviews, vol. 15,
% no. 3 (August), pp. 197-235, revised version of Brock et al. (1987)
%
% Dechert, Davis (1988), "BDS STATS: A Program to Calculate the Statistics of the
% Grassberger-Procaccia Correlation Dimension Based on the Paper "A Test for
% Independence" by W. A. Brock, W. D. Dechert and J. A. Scheinkman", version 8.21
% (latest), MS-DOS software available on gopher.econ.wisc.edu
%
% Kanzler, Ludwig (1998), "Very Fast and Correctly Sized Estimation of the BDS Statistic",
% Oxford University, Department of Economics, working paper, available on
% http://users.ox.ac.uk/~econlrk
%
% LeBaron, Blake (1988, 1990, 1997a), "BDSTEST.C", version June 1997 (latest), C source
% code available on gopher.econ.wisc.edu
%
% LeBaron, Blake (1997b), "A Fast Algorithm for the BDS Statistic", Studies in
% Nonlinear Dynamics and Econometrics, vol. 2, pp. 53-59
% ACKNOWLEDGEMENT:
%
% I am grateful to Blake LeBaron for giving me the exclusive opportunity to beta-test his
% C programme in its compiled version for MATLAB 5 and thus enabling me to compare the two
% programmes directly. I have benefited from the many associated discussions.
% End of file.
end