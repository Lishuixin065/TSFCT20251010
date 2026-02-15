%% Algorithm for multiple compressions - Optimized Version
function [tfrsq1] = Msqueeze_FCT(x, tfc0, R0, T0, t, tcrtic, num, thre)
%
% Optimized version with vectorization and pre-computation
%
% Inputs:
%   x       : Input signal (vector, Nx1)
%   tfc0    : Initial Frequency-domain Chirplet Transform (FCT) 
%   R0      : Group Delay Dispersion (GDD) reassignment operator
%   T0      : Group Delay (GD) reassignment operator 
%   t       : Time axis vector (1xN)
%   tcrtic  : GDD axis vector 
%   num     : Number of squeezing iterations (positive integer)
%   thre    : Energy concentration threshold (0 < thre < 1)
%
% Output:
%   tfrsq1  : HTSFCT
%

[a, b, c] = size(tfc0); 
Threshold = thre * mean(abs(x).^2);

% 预计算步长
del_t = t(2) - t(1);
del_c = tcrtic(2) - tcrtic(1);

% 预计算R和T的索引 - 使用向量化操作


% 向量化计算R和T（比双循环快）
R = round(real(R0) ./ del_c) + floor(length(tcrtic)/2) + 1;
T = round(real(T0) ./ del_t) + 1;

% 边界检查
R(R < 1) = 1;
R(R > a) = a;
T(T < 1) = 1;
T(T > b) = b;

fprintf('Total frequency bins: %d; Current:     ', c);

tfrsq1 = zeros(a, b, c);

% 主循环 - 提前提取所有频率切片
for kk = 1:num
    tfrsq_current = zeros(a, b, c);
    
    for fidx = 1:c
        fprintf('\b\b\b\b');
        fprintf('%4d', fidx);
        
        % 提前获取所有需要的切片
        tf0_slice = tfc0(:, :, fidx);
        R_slice = R(:, :, fidx);
        T_slice = T(:, :, fidx);
        
        % 创建稀疏累加矩阵 - 使用accumarray进行向量化累加
        sst = zeros(a, b);
        
        % 找出所有超过阈值的点
        valid_mask = abs(tf0_slice) > Threshold;
        
        if any(valid_mask(:))
            % 获取有效点的索引和值
            [valid_rows, valid_cols] = find(valid_mask);
            valid_values = tf0_slice(valid_mask);
            
            % 获取对应的目标位置
            dest_rows = R_slice(valid_mask);
            dest_cols = T_slice(valid_mask);
            
            % 使用线性索引进行累加（比双重循环快）
            % 将二维目标位置转换为线性索引
            dest_linear = sub2ind([a, b], dest_rows, dest_cols);
            
            % 使用accumarray进行累加（更高效）
            sst_linear = accumarray(dest_linear, valid_values, [a*b, 1]);
            sst = reshape(sst_linear, [a, b]);
        end
        
        tfrsq_current(:, :, fidx) = sst;
    end
    
    % 更新tfc0用于下一次迭代
    if kk < num
        tfc0 = tfrsq_current;
    end
    
    tfrsq1 = tfrsq_current;
end

fprintf('\n');
end