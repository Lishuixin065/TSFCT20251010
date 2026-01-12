function [tfrsq1] = Msqueeze_CT_fixed(x, ct, lambda, omega, tfrtic, tcrtic, num)
    % 修复版压缩变换 - 保持与原代码一致的逻辑
    
    % 获取尺寸
    [a, b, c] = size(ct);
    
    % 预计算常数
    Threshold = 0.0001 * mean(abs(x(:)).^2);
    del_c = tcrtic(2) - tcrtic(1);
    del_f = tfrtic(2) - tfrtic(1);
    
    fprintf('预处理阶段...\n');
    
    % 转换为单精度减少内存使用
    ct = single(ct);
    lambda = single(lambda);
    omega = single(omega);
    
    fprintf('开始压缩处理...\n');
    
    % 预分配输出
    tfrsq1 = zeros(a, b, c, 'single');
    
    % 进度显示初始化
    fprintf('总时间迭代: %d; 当前:     ', c);
    
    for kk = 1:num
        for tidx = 1:c
            % 进度显示
            if mod(tidx, max(1, round(c/20))) == 0 || tidx == c
                fprintf('\b\b\b\b\b\b\b\b\b');
                tmp = sprintf('%4d/%4d', tidx, c);
                fprintf(tmp);
            end
            
            % 提取当前时间片
            lambda2 = squeeze(lambda(:, :, tidx));
            omega2 = squeeze(omega(:, :, tidx));
            tf0 = squeeze(ct(:, :, tidx));
            
            % 归一化坐标
            omega2 = round(omega2 / del_f) + 1;
            lambda2 = round((lambda2 + tcrtic(end)) / del_c) + 1;
            
            % 阈值处理
            invalid_mask = (abs(tf0) < Threshold) | (abs(omega2) < Threshold) | (abs(lambda2) < Threshold);
            omega2(invalid_mask) = 0;
            lambda2(invalid_mask) = 0;
            
            % 创建有效掩码
            valid_mask = (abs(tf0) > Threshold) & ...
                        (omega2 >= 1) & (omega2 <= b) & ...
                        (lambda2 >= 1) & (lambda2 <= a);
            
            if any(valid_mask(:))
                % 提取有效数据
                [row_idx, col_idx] = find(valid_mask);
                
                if ~isempty(row_idx)
                    % 提取对应的值
                    k_vals = omega2(valid_mask);
                    m_vals = lambda2(valid_mask);
                    tf_vals = tf0(valid_mask);
                    
                    % 过滤有效范围内的索引
                    valid_range = (k_vals >= 1) & (k_vals <= b) & ...
                                 (m_vals >= 1) & (m_vals <= a);
                    
                    if any(valid_range)
                        k_vals = k_vals(valid_range);
                        m_vals = m_vals(valid_range);
                        tf_vals = tf_vals(valid_range);
                        
                        % 保持与原代码一致的稀疏矩阵操作
                        % 转换为double类型用于sparse函数
                        sst_sparse = sparse(double(m_vals), double(k_vals), double(tf_vals), a, b);
                        sst = single(full(sst_sparse));
                    else
                        sst = zeros(a, b, 'single');
                    end
                else
                    sst = zeros(a, b, 'single');
                end
            else
                sst = zeros(a, b, 'single');
            end
            
            tfrsq1(:, :, tidx) = sst;
        end
        
        % 更新用于下一次迭代
        if kk < num
            ct = tfrsq1;
        end
    end
    
    fprintf('\n压缩完成!\n');
end



