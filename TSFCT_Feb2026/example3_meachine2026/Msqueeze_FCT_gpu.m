function [tfc0,tfrsq1,tfrtic,aa,R,T] = Msqueeze_FCT_gpu(x,Hz,sigma,arange,hop,num,thre)

[~,xcol] = size(x);
if (xcol~=1)
    error('x must have only one column');
end

N = length(x);
t0 = 0:hop:N-1;
t = t0/Hz;
del_t = hop/Hz;

if mod(N,2)==0
    tfrtic = (0:hop:N/2)*Hz/N;
else
    tfrtic = (0:hop:(N-1)/2)*Hz/N;
end

f = tfrtic;

tLen = length(t);
fLen = length(tfrtic);
ht = -tLen:tLen;
Lh = tLen;
time = ht/Hz;

aa = linspace(-arange, arange, 2*round(tLen/2)+1);
aLen = length(aa);

% GPU检查
if gpuDeviceCount == 0
    error('No GPU available for this function');
end

gpuDevice(1);

% 所有数据转移到GPU
x_gpu = gpuArray(single(x));
aa_gpu = gpuArray(single(aa));
time_gpu = gpuArray(single(time));
f_gpu = gpuArray(single(f));
t_gpu = gpuArray(single(t));
t0_gpu = gpuArray(single(t0));

% 预分配GPU结果数组
tfc0_gpu = gpuArray.zeros(aLen, tLen, fLen, 'single');
tfc1_gpu = gpuArray.zeros(aLen, tLen, fLen, 'single');
tfc2_gpu = gpuArray.zeros(aLen, tLen, fLen, 'single');
T_gpu = gpuArray.zeros(aLen, tLen, fLen, 'single');
R_gpu = gpuArray.zeros(aLen, tLen, fLen, 'single');

fprintf('GPU vectorized computation. Time frames total: %d; now:     ', tLen);

% 完全保持原代码逻辑，只是在GPU上执行
for tidx = 1:tLen
    if mod(tidx, 10) == 0
        fprintf('\b\b\b\b');
        fprintf('%4d', tidx);
    end
    
    ti = t0_gpu(tidx) + 1;
    
    % 计算tau范围和索引（在GPU上计算）
    tau_min = -min([round(tLen/2)-1, tLen, ti-1]);
    tau_max = min([round(tLen/2)-1, tLen, N-ti]);
    tau = tau_min:tau_max;
    
    if isempty(tau)
        continue;
    end
    
    indices = rem(tLen + tau, tLen) + 1;
    tau_len = length(tau);
    
    % 获取信号片段
    x_segment = x_gpu(ti + tau).';
    
    % 向量化处理所有chirp率
    L0 = single(1) + 1i * single(2 * pi) * aa_gpu(:) .* single(sigma^2);
    
    % 计算窗口函数值
    time_subset = time_gpu(Lh + 1 - tau);
    
    % 扩展维度
    L0_expanded = repmat(L0, 1, tau_len);
    time_expanded = repmat(time_subset, aLen, 1);
    
    % 计算窗口函数
    exp_factor = exp(-single(2*pi^2) * single(sigma^2) * time_expanded.^2 ./ L0_expanded);
    
    h = (single(1) ./ sqrt(L0_expanded)) .* exp_factor;
    th = (-1i * single(2 * pi) * single(sigma^2) * time_expanded .* (L0_expanded.^(-single(3/2)))) .* exp_factor;
    tth = single(sigma^2) .* ((L0_expanded.^(-single(3/2))) - ((single(2*pi*sigma) * time_expanded).^2) .* (L0_expanded.^(-single(5/2)))) .* exp_factor;
    
    % 共轭
    h = conj(h);
    th = conj(th);
    tth = conj(tth);
    
    % 信号扩展
    x_expanded = repmat(x_segment, aLen, 1);
    
    % 计算TF segments
    tf0_seg = x_expanded .* h;
    tf1_seg = x_expanded .* th;
    tf2_seg = x_expanded .* tth;
    
    % 构建完整TF矩阵
    tf0_full = complex(gpuArray.zeros(aLen, tLen, 'single'));
    tf1_full = complex(gpuArray.zeros(aLen, tLen, 'single'));
    tf2_full = complex(gpuArray.zeros(aLen, tLen, 'single'));
    
    % 填充数据
    tf0_full(:, indices) = tf0_seg;
    tf1_full(:, indices) = tf1_seg;
    tf2_full(:, indices) = tf2_seg;
    
    % FFT
    tf0 = fft(tf0_full, [], 2);
    tf1 = fft(tf1_full, [], 2);
    tf2 = fft(tf2_full, [], 2);
    
    % 只取正频率
    tf0 = tf0(:, 1:fLen);
    tf1 = tf1(:, 1:fLen);
    tf2 = tf2(:, 1:fLen);
    
    % 相位校正
    phase_corr = exp(-single(2*pi) * 1i * f_gpu(:).' * t_gpu(tidx));
    tf0 = tf0 .* repmat(phase_corr, aLen, 1);
    tf1 = tf1 .* repmat(phase_corr, aLen, 1);
    tf2 = tf2 .* repmat(phase_corr, aLen, 1);
    
    % 存储结果
    tfc0_gpu(:, tidx, :) = tf0;
    tfc1_gpu(:, tidx, :) = tf1;
    tfc2_gpu(:, tidx, :) = tf2;
    
    % 重分配计算
    M0 = tf2 .* tf0 - tf1 .* tf1;
    M1 = -tf0 .* tf1;
    M2 = tf0 .* tf0;
    
    % 阈值计算
    avg_M0 = mean(abs(M0), 2);
    thre = single(0.001) * avg_M0;
    E0 = abs(M0) > repmat(thre, 1, fLen);
    
    % 计算重分配参数
    omega1 = single(tidx-1) * del_t + M1 ./ (M0 + (M0==0)) ./ (single(2.0*pi) * 1i);
    lambda1 = repmat(aa_gpu(:), 1, fLen) + M2 ./ (M0 + (M0==0)) ./ (single(2.0*pi) * 1i);
    
    % 应用阈值
    R_gpu(:, tidx, :) = lambda1 .* single(E0);
    T_gpu(:, tidx, :) = omega1 .* single(E0);
end


% 将tcrtic也上传到GPU
tcrtic_gpu = gpuArray(single(aa));

% 调用Msqueeze_FCT（假设它已支持GPU或使用gpuArray输入）
[tfrsq1_gpu] =  Msqueeze_FCT_gpu_vectorized(x_gpu, tfc0_gpu, R_gpu, T_gpu, t_gpu, tcrtic_gpu, num, thre);

% 将结果从GPU传回CPU
tfc0 = gather(tfc0_gpu);
tfrsq1 = gather(tfrsq1_gpu);

T = gather(T_gpu);
R = gather(R_gpu);

fprintf('GPU computation completed!\n');
end

function [tfrsq1] = Msqueeze_FCT_gpu_vectorized(x, tfc0, R0, T0, t, tcrtic, num, thre)
% GPU向量化版本

[a, b, c] = size(tfc0); 
Threshold = thre * mean(abs(x).^2);

% 预计算步长
del_t = t(2) - t(1);
del_c = tcrtic(2) - tcrtic(1);

% GPU计算
R = round(real(R0) ./ del_c) + floor(length(tcrtic)/2) + 1;
T = round(real(T0) ./ del_t) + 1;

% 边界检查
R = max(1, min(a, R));
T = max(1, min(b, T));

fprintf('Total frequency bins: %d; Processing...\n', c);

% 初始化结果
tfrsq1 = gpuArray.zeros(a, b, c, 'single');

for kk = 1:num
    tfrsq_current = gpuArray.zeros(a, b, c, 'single');
    
    % 对于每个频率切片
    for fidx = 1:c
        if mod(fidx, 100) == 0
            fprintf('Processing slice %d/%d\n', fidx, c);
        end
        
        tf0_slice = tfc0(:, :, fidx);
        R_slice = R(:, :, fidx);
        T_slice = T(:, :, fidx);
        
        % 创建有效点索引
        valid_mask = abs(tf0_slice) > Threshold;
        num_valid = nnz(valid_mask);
        
        if num_valid > 0
            % 提取有效数据
            valid_values = tf0_slice(valid_mask);
            
            % 获取目标位置
            dest_positions = [R_slice(valid_mask), T_slice(valid_mask)];
            
            % 确保目标位置在范围内
            dest_positions(:, 1) = max(1, min(a, dest_positions(:, 1)));
            dest_positions(:, 2) = max(1, min(b, dest_positions(:, 2)));
            
            % 累加到结果
            sst = gpuArray.zeros(a, b, 'single');
            
            % GPU累加：使用线性索引
            linear_indices = sub2ind([a, b], dest_positions(:, 1), dest_positions(:, 2));
            
            % 使用unique和accumarray的替代方案
            [unique_indices, ~, ic] = unique(linear_indices);
            accumulated_values = accumarray(ic, double(valid_values));
            
            % 将累加结果放回矩阵
            sst(unique_indices) = accumulated_values;
            
            tfrsq_current(:, :, fidx) = sst;
        end
    end
    
    if kk < num
        tfc0 = tfrsq_current;
    end
    
    tfrsq1 = tfrsq_current;
end

fprintf('GPU squeezing completed!\n');
end