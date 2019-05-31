function [ output_args ] = computeCapacity( a, b )
%计算CUE用户的容量，根据文章中的定理2
%a/((a-b)*log2)*(exp(1/a)*expint(1/a) - exp(1/b)*expint(1/b))


if a>=(1/700) && b>=(1/700)
    output_args = a/((a-b)*log(2))*(exp(1/a)*expint(1/a) - exp(1/b)*expint(1/b));%expint：指数积分。
elseif a<(1/700) && b<(1/700)
    output_args = a/((a-b)*log(2))*(a - b);
elseif b < (1/700)
    output_args = a/((a-b)*log(2))*(exp(1/a)*expint(1/a) - b); % 发现：exp(x)*expint(x)= 1/x当x无穷大时成立
elseif a < (1/700)
    output_args = a/((a-b)*log(2))*(a - exp(1/b)*expint(1/b));
end

end

