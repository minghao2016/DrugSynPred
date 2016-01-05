%binary_superset binary_subset
%using decimal value to represent the super set and subset
% a = [1,0,0,1];
% v1 superset does not contain itself
function [b_superset,b_subset]=binary_set(bs)
lbs=length(bs);
dec_zeros = 2.^(lbs-find(bs==0));
dec_ones = 2.^(lbs-find(bs==1));
dec_bs = 2^lbs-1-sum(dec_zeros);


zeros_gcode = dec2bin(grays_2(length(dec_zeros)))-'0';
% true superset remove gray_code first line
% zeros_gcode(1,:)=[];

ones_gcode = dec2bin(grays_2(length(dec_ones)))-'0';
ones_gcode(1,:)=[];

dec_zeros_bs = repmat(dec_zeros,length(zeros_gcode),1);

if(isempty(dec_zeros_bs))
    b_superset=2^lbs-1;
    else
    super_gcode = dec_zeros_bs.*zeros_gcode;
    superset_sum = sum(super_gcode,2);
    b_superset = superset_sum(:)+ dec_bs;
end

dec_ones_bs = repmat(dec_ones,length(ones_gcode),1);

if(isempty(dec_ones_bs))
    b_subset=[];
else
    subset_sum =  sum( dec_ones_bs.*ones_gcode,2);
    b_subset = dec_bs - subset_sum(:) ;
end

function  G = grays_2(n)
%  G = grays(n)  is a column of  2^n  distinct  n-bit
%  integers that step through the  Gray Codes  running
%  from  G(1) = 000...000  to  G(2^n) = 100...000 .
%  Each  G(j+1)  is obtained from  G(j)  by changing
%  just one bit of  G(j) .  To see which bit,  compute
%  whichbit = bitxor(G(j), G(j+1))  to get  2^m  for a
%  nonnegative integer  m < n .  Keep  n < 27  because
%  grays(n)  costs time and memory proportional to  2^n .
%  See also  graystep.m,  int2gray.m,  gray2int.m .
%  Display  G + 2^52  in  hex  to see how its last  n
%  bits change.                 W. Kahan,  8 July 2007

n = n(:) ;  T = length(n) ;
% if ( (T~=1)|(n~=round(n))|(n<1)|(n>26) ),  N = n ,
%   error(' grays(N)  needs a small positive integer  N .'),  
% end
G = zeros(2^n,1) ;  G(2) = 1 ;  T = 2 ;
for  k = 2:n
    T2 = T+T ;  ... = 2^k
    G(T+1:T2) = T + flipud(G(1:T)) ;
    T = T2 ;  
end
