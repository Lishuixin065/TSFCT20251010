function [RQF_IF, RQF_CR] = RQF_ridge(if1, if2, chirp1, chirp2, der_x1, dder_x1, der_x2, dder_x2)
% RQF_RIDGE 计算瞬时频率和chirp率的重构质量因子RQF
% 
% 输入参数:
%   if1, if2     : 估计的瞬时频率 (两个分量)
%   chirp1, chirp2 : 估计的chirp率 (两个分量)
%   der_x1, der_x2 : 真实的瞬时频率 (信号x1和x2)
%   dder_x1, dder_x2 : 真实的chirp率 (信号x1和x2)
% 
% 输出参数:
%   RQF_IF1 : 信号x1瞬时频率的RQF [dB]
%   RQF_IF2 : 信号x2瞬时频率的RQF [dB]
%   RQF_CR1 : 信号x1 chirp率的RQF [dB]
%   RQF_CR2 : 信号x2 chirp率的RQF [dB]
% 
% 公式: RQF = 10*log10(||真值||² / ||真值 - 估计值||²)

% 1. 边界处理（与你原代码完全一致）
len = length(der_x1);
tt = round(len/8) : round(7*len/8);


% 2. 分量匹配与RQF计算
% 瞬时频率部分
if sum(abs(if1(tt) - der_x1(tt)).^2) < sum(abs(if2(tt) - der_x1(tt)).^2)
    % if1更接近der_x1, if2更接近der_x2
    IF_x1_true = der_x1(tt);
    IF_x1_est = if1(tt);
    IF_x2_true = der_x2(tt);
    IF_x2_est = if2(tt);
else
    % if2更接近der_x1, if1更接近der_x2
    IF_x1_true = der_x1(tt);
    IF_x1_est = if2(tt);
    IF_x2_true = der_x2(tt);
    IF_x2_est = if1(tt);
end

% 计算瞬时频率的RQF
norm_IF1_true = norm(IF_x1_true, 2)^2;
norm_IF1_err = norm(IF_x1_true - IF_x1_est, 2)^2;
RQF_IF1 = 10 * log10(norm_IF1_true / norm_IF1_err);

norm_IF2_true = norm(IF_x2_true, 2)^2;
norm_IF2_err = norm(IF_x2_true - IF_x2_est, 2)^2;
RQF_IF2 = 10 * log10(norm_IF2_true / norm_IF2_err);

% chirp率部分
if sum(abs(chirp1(tt) - dder_x1(tt)).^2) < sum(abs(chirp2(tt) - dder_x1(tt)).^2)
    % chirp1更接近dder_x1, chirp2更接近dder_x2
    CR_x1_true = dder_x1(tt);
    CR_x1_est = chirp1(tt);
    CR_x2_true = dder_x2(tt);
    CR_x2_est = chirp2(tt);
else
    % chirp2更接近dder_x1, chirp1更接近dder_x2
    CR_x1_true = dder_x1(tt);
    CR_x1_est = chirp2(tt);
    CR_x2_true = dder_x2(tt);
    CR_x2_est = chirp1(tt);
end

% 计算chirp率的RQF
norm_CR1_true = norm(CR_x1_true, 2)^2;
norm_CR1_err = norm(CR_x1_true - CR_x1_est, 2)^2;
RQF_CR1 = 10 * log10(norm_CR1_true / norm_CR1_err);

norm_CR2_true = norm(CR_x2_true, 2)^2;
norm_CR2_err = norm(CR_x2_true - CR_x2_est, 2)^2;
RQF_CR2 = 10 * log10(norm_CR2_true / norm_CR2_err);

RQF_IF=RQF_IF1+RQF_IF2;
RQF_CR=RQF_CR1+RQF_CR2;

end