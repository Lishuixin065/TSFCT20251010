function [RQF] = RQF_recov(x1, x2, recov)
% RQF_RECOV_ERROR 计算两个信号的重构质量因子RQF
% 
% 输入参数:
%   x1   : 原始信号1 (1×N行向量)
%   x2   : 原始信号2 (1×N行向量)
%   recov: 重构信号 (2×N矩阵，两行分别对应两个重构信号)
% 
% 输出参数:
%   RQF1 : 信号1的RQF值 [dB]
%   RQF2 : 信号2的RQF值 [dB]
% 
% 公式: RQF = 10*log10(||x||² / ||x - x_recov||²)

% 1. 边界处理（与你原代码完全一致）
tn = length(x1);
time = ceil(tn/8) : floor(tn*7/8);


% 提取核心时段
x1_core = x1(1, :);
x2_core = x2(1, :);
recov_core = recov(:,:);

% 2. 分量匹配判断（与你原代码逻辑相同）
% 判断哪个重构分量对应哪个原始信号
if sum(abs(x1_core - recov_core(1, :)), 'all') < sum(abs(x2_core - recov_core(1, :)), 'all')
    % 情况1: recov(1,:) 更接近 x1, recov(2,:) 更接近 x2
    x1_recov = recov_core(1, :);
    x2_recov = recov_core(2, :);
else
    % 情况2: recov(1,:) 更接近 x2, recov(2,:) 更接近 x1
    x1_recov = recov_core(2, :);
    x2_recov = recov_core(1, :);
end

% 3. 计算RQF（使用你提供的公式）
% 信号1的RQF
norm_x1_sq = norm(x1_core, 2)^2;
norm_err1_sq = norm(x1_core - x1_recov, 2)^2;
RQF1 = 10 * log10(norm_x1_sq / norm_err1_sq);

% 信号2的RQF
norm_x2_sq = norm(x2_core, 2)^2;
norm_err2_sq = norm(x2_core - x2_recov, 2)^2;
RQF2 = 10 * log10(norm_x2_sq / norm_err2_sq);
RQF=RQF1+RQF2;
end