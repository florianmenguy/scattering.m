function y_ft = multiply_fft(x_ft,filter_struct,log2_resampling,colons,subscripts)
%% In-place Fourier transform
for subscript = subscripts
    x_ft = fft(x_ft,[],subscript);
end

%% Definition of x_sizes
if length(subscripts)>1
    error('Fourier-based filtering not ready in dimension >1');
end
input_sizes = size(x_ft);
x_sizes = input_sizes(subscripts);

%% Loading of filter
filter_ft = filter_struct.ft;
filter_start = filter_struct.ft_start;

%% Trimming of x_ft if needed
y_sizes = pow2(x_sizes,log2_resampling);
pos_range_start = max(1,filter_start);
x_end = x_sizes/2;
filter_sizes = length(filter_ft);
filter_end =  filter_start + filter_sizes - 1;
filter_mod_end = mod(filter_end-1,x_sizes) + 1;
y_end = y_sizes/2;
pos_range_end = min([x_end,filter_mod_end,y_end]);
x_start = x_sizes/2 + 1;
filter_mod_start = mod(filter_start-1,x_sizes) + 1;
y_start = x_sizes - y_sizes/2 + 1;
neg_range_start = max([x_start,filter_mod_start,y_start]);
if filter_mod_start>x_end
    neg_range_end = x_sizes;
elseif filter_mod_end>x_end
    neg_range_end = min([x_sizes,filter_mod_end]);
else
    neg_range_end = 1;
end
pos_x_range = pos_range_start:pos_range_end;
neg_x_range = neg_range_start:neg_range_end;
x_range = cat(2,neg_x_range,pos_x_range);
colons.subs{subscripts} = x_range;
trimmed_x_ft = subsref(x_ft,colons);

%% Trimming of filter_ft if needed
pos_filter_range = mod(pos_x_range - filter_start,filter_sizes) + 1;
unbounded_neg_filter_range = ...
    neg_x_range - filter_start + max(filter_sizes-x_sizes,0);
neg_filter_range = mod(unbounded_neg_filter_range,filter_sizes) + 1;
filter_range = cat(2,neg_filter_range,pos_filter_range);
colons.subs{subscripts} = filter_range;
trimmed_filter_ft = subsref(filter_ft,colons);

%% Product between non-negligible coefficients of x_ft and filter_ft
neg_y_range_start = neg_range_start + y_sizes - x_sizes;
neg_y_range_end = neg_range_end + y_sizes - x_sizes;
y_range = cat(2,neg_y_range_start:neg_y_range_end,pos_x_range);
colons.subs{subscripts} = y_range;
sub_y_ft = bsxfun(@times,trimmed_x_ft,trimmed_filter_ft);

%% Assignment of product to zero-allocated y_ft
y_tensor_sizes = size(x_ft);
y_tensor_sizes(subscripts) = y_sizes;
y_ft = zeros(y_tensor_sizes);
y_ft = subsasgn(y_ft,colons,sub_y_ft);
end