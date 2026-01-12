function [tfrsq1] = Msqueeze0_optimized(x, tfc0, R0, T0, t, tcrtic, num, thre)
% 优化版本 - 使用向量化、预计算和内存优化
%
% 主要优化点：
% 1. 预计算所有索引
% 2. 使用accumarray的向量化版本
% 3. 减少内存分配和拷贝
% 4. 并行化处理

[a, b, c] = size(tfc0); 
Threshold = thre * mean(abs(x).^2);

% 预计算步长
del_t = t(2) - t(1);
del_c = tcrtic(2) - tcrtic(1);

% ==================== 优化1：向量化索引计算 ====================
fprintf('预计算索引...\n');
tic;

% 一次性计算所有索引（比在循环中计算快）
R = min(max(round(real(R0) ./ del_c) + floor(length(tcrtic)/2) + 1, 1), a);
T = min(max(round(real(T0) ./ del_t) + 1, 1), b);

% 转换为线性索引（加速后续访问）
R_linear = reshape(R, a*b*c, 1);
T_linear = reshape(T, a*b*c, 1);

% 预计算目标线性索引
dest_linear_all = sub2ind([a, b], R_linear, T_linear);
dest_linear_all = reshape(dest_linear_all, a, b, c);

fprintf('索引预计算完成，耗时: %.3f秒\n', toc);

% ==================== 优化2：预计算有效点掩码 ====================
fprintf('预计算阈值掩码...\n');
tic;

% 一次性计算所有有效点
tfc0_linear = reshape(tfc0, a*b*c, 1);
valid_mask_linear = abs(tfc0_linear) > Threshold;

% 预分配工作空间
tfrsq1 = zeros(a, b, c);



fprintf('开始压缩迭代...\n');


% ==================== 优化3：主循环优化 ====================
for kk = 1:num
    fprintf('迭代 %d/%d: ', kk, num);
    
    % 预分配当前迭代结果
    tfrsq_current = zeros(a, b, c);
    
    % ==================== 优化4：分块并行处理 ====================
    % 每个频率切片独立处理，适合并行化
    
    % 如果没有并行工具箱，使用普通循环
    % 如果有并行工具箱，可以改为 parfor
    
    % 进度显示
    progress_interval = max(1, floor(c/20));
    
    for fidx = 1:c
        % 进度显示
        if mod(fidx, progress_interval) == 0
            fprintf('.');
        end
        
        % 获取当前切片的线性索引
        slice_start = (fidx-1)*a*b + 1;
        slice_end = fidx*a*b;
        
        % 提取当前切片的数据
        slice_valid_mask = valid_mask_linear(slice_start:slice_end);
        
        if any(slice_valid_mask)
            % 提取有效点的值和目标位置
            valid_values = tfc0_linear(slice_start:slice_end);
            valid_values = valid_values(slice_valid_mask);
            
            dest_linear = dest_linear_all(:, :, fidx);
            dest_linear = dest_linear(:);
            dest_linear = dest_linear(slice_valid_mask);
            
            % ==================== 优化5：快速累加 ====================
            % 方法1：使用accumarray（标准方法）
            sst_linear = accumarray(dest_linear, valid_values, [a*b, 1]);
            
            % 方法2：对于稀疏数据，使用sparse更快
            % if nnz(slice_valid_mask) < a*b/100  % 少于1%的点
            %     sst_linear = full(sparse(dest_linear, 1, valid_values, a*b, 1));
            % else
            %     sst_linear = accumarray(dest_linear, valid_values, [a*b, 1]);
            % end
            
            % 重塑为二维矩阵
            tfrsq_current(:, :, fidx) = reshape(sst_linear, [a, b]);
        end
    end
    
 
    % ==================== 优化6：更新数据用于下一次迭代 ====================
    if kk < num
        tfc0 = tfrsq_current;
        tfc0_linear = reshape(tfc0, a*b*c, 1);
        
        % 重新计算有效点掩码（阈值可能变化）
        valid_mask_linear = abs(tfc0_linear) > Threshold;
    end
    
    tfrsq1 = tfrsq_current;
end

end