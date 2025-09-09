%****************************************************************************************
%*                                                                                      *
%*                      EE 388 - Spatial Coupling LDPC (SNR Version)                   *
%*                                                                                      *
%****************************************************************************************

function spatial_coupling_ldpc_simulation_snr()
    % 参数定义
    m_PROTOGRAPH = 100;
    L_SPATIALCOUPLING = 500;
    k_LDPC_reg_ensemble = 10;
    l_LDPC_reg_ensemble = 5;    % k/l = q, 必须是整数
    
    MAX_BP_ITER = 200;
    MAX_BP_UNCHANGED_ITER = 10;
    SIMULATION_NUM = 20;
    
    % SNR范围 (dB)
    snr_vals = [-2, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    
    fprintf('开始空间耦合LDPC SNR仿真...\n');
    
    % 代码构造
    m = m_PROTOGRAPH;
    l = l_LDPC_reg_ensemble;
    k = k_LDPC_reg_ensemble;
    L = L_SPATIALCOUPLING;
    q = k / l;
    
    fprintf('参数: m=%d, l=%d, k=%d, L=%d, q=%d\n', m, l, k, L, q);
    
    varNum = q * L * m;
    chkNum = chkSize_SC_regLDPC(l, m, L);
    
    fprintf('变量节点数: %d, 校验节点数: %d\n', varNum, chkNum);
    
    % 初始化连接矩阵
    varConnectOriginal = create_connection_matrix(q*L, l+1);
    
    % 构建耦合结构
    fprintf('构建耦合结构...\n');
    varConnectOriginal = coupling_regular(l, q, L, varConnectOriginal);
     
    % 构建原型图
    fprintf('构建原型图...\n');
    varConnect = create_connection_matrix(varNum, l+1);
    varConnect = protograph(m, chkNum/m, q*L, varConnectOriginal, varConnect);
    
    % 完成图结构
    fprintf('完成图结构...\n');
    [varConnect, chkConnect, degVar, degChk, varConnectPosition, chkConnectPosition] = ...
        complete_graph(varNum, chkNum, k, l, varConnect);
    
    fprintf('代码构造完成.\n');
    
    % 仿真参数初始化
    BER = zeros(1, length(snr_vals));
    FER = zeros(1, length(snr_vals));  % 帧错误率
    
    progress = 0;
    total_sims = SIMULATION_NUM * length(snr_vals);
    
    % 主仿真循环
    parfor cnt_snr = 1:length(snr_vals)
        snr_db = snr_vals(cnt_snr);
        snr_linear = 10^(snr_db/10);  % 转换为线性值
        
        ber_sum = 0;
        fer_sum = 0;
        
        for cnt_simul = 1:SIMULATION_NUM
            % 生成随机信息位
            info_bits = randi([0, 1], varNum, 1);
            
            % 信道仿真 - AWGN信道
            [received_signal, llr] = AWGN_channel(info_bits, snr_linear);
            
            % BP解码
            decodedword = BPDecoder_AWGN(...
                llr, varConnect, chkConnect, varConnectPosition, chkConnectPosition, ...
                degVar, degChk, varNum, chkNum, l+1, k+1, MAX_BP_ITER, MAX_BP_UNCHANGED_ITER);
            
            % 计算BER
            error_bits = sum(decodedword ~= info_bits);
            ber_sum = ber_sum + error_bits / varNum;
            
            % 计算FER (如果有一个或多个错误位，则认为帧错误)
            if error_bits > 0
                fer_sum = fer_sum + 1;
            end
        end
        
        BER(cnt_snr) = ber_sum / SIMULATION_NUM;
        FER(cnt_snr) = fer_sum / SIMULATION_NUM;
        fprintf('SNR = %.1f dB, BER = %.10f, FER = %.10f\n', snr_db, BER(cnt_snr), FER(cnt_snr));
    end
    
    % 显示最终结果
    fprintf('\n最终结果:\n');
    for cnt_snr = 1:length(snr_vals)
        snr_db = snr_vals(cnt_snr);
        fprintf('SNR = %.1f dB, BER = %.10f, FER = %.10f\n', snr_db, BER(cnt_snr), FER(cnt_snr));
    end
    
    % 保存结果到文件
    filename = sprintf('Result_SNR_%d.txt', L_SPATIALCOUPLING);
    save_results_snr(filename, snr_vals, BER, FER);
    
    % 绘图
    plot_results_snr(snr_vals, BER, FER);
    
    fprintf('SNR仿真完成!\n');
end

%% 辅助函数

% 创建Fisher-Yates随机排列
function vector = randperm_custom(n)
    vector = 0:(n-1);  % MATLAB使用1-based索引，但这里保持0-based逻辑
    for k = 1:n
        j = randi([k, n]);
        temp = vector(j);
        vector(j) = vector(k);
        vector(k) = temp;
    end
end

% 创建连接矩阵
function matrix = create_connection_matrix(rows, cols)
    matrix = -ones(rows, cols, 'int32');
end

% 创建L耦合正则LDPC码
function varConnect = coupling_regular(l, q, L, varConnect)
    if mod(l, 2) == 1
        l1 = (l - 1) / 2;
        l2 = (l - 1) / 2;
    else
        l1 = (l - 2) / 2;
        l2 = l / 2;
    end
    
    for i = 1:L
        cnt = 1;  % MATLAB 1-based索引
        for a = i:(i + l - 1)
            for r = 1:q
                row_idx = (r-1)*L + i;  % 转换为MATLAB索引
                if row_idx <= size(varConnect, 1) && cnt <= size(varConnect, 2)
                    varConnect(row_idx, cnt) = a - 1;  % 保持0-based的值
                end
            end
            cnt = cnt + 1;
        end
    end
end

% 创建m倍原型图代码
function varConnect = protograph(m, p, n, varConnectProtograph, varConnect)
    vector = zeros(1, m);
    posVar = zeros(1, m*n);
    
    for i = 1:n
        cnt = 1;
        while cnt <= size(varConnectProtograph, 2)
            if varConnectProtograph(i, cnt) ~= -1
                a = varConnectProtograph(i, cnt);
                vector = randperm_custom(m);
                for b = 1:m
                    row_idx = (i-1)*m + b;
                    col_idx = posVar(row_idx) + 1;
                    if row_idx <= size(varConnect, 1) && col_idx <= size(varConnect, 2)
                        varConnect(row_idx, col_idx) = vector(b) + a*m;
                        posVar(row_idx) = posVar(row_idx) + 1;
                    end
                end
                cnt = cnt + 1;
            else
                break;
            end
        end
    end
end

% 完成图结构并去除重复边
function [varConnect, chkConnect, degVar, degChk, varConnectPosition, chkConnectPosition] = ...
    complete_graph(varNum, chkNum, k, l, varConnect)
    
    % 初始化输出矩阵
    chkConnect = create_connection_matrix(chkNum, k+1);
    varConnectPosition = create_connection_matrix(varNum, l+1);
    chkConnectPosition = create_connection_matrix(chkNum, k+1);
    degVar = zeros(varNum, 1);
    degChk = zeros(chkNum, 1);
    
    % 去除重复边
    for i = 1:varNum
        unique_connections = [];
        for cnt = 1:size(varConnect, 2)
            if varConnect(i, cnt) == -1
                break;
            end
            if ~ismember(varConnect(i, cnt), unique_connections)
                unique_connections = [unique_connections, varConnect(i, cnt)];
            end
        end
        
        % 更新连接
        for j = 1:length(unique_connections)
            if j <= size(varConnect, 2)
                varConnect(i, j) = unique_connections(j);
            end
        end
        for j = (length(unique_connections)+1):size(varConnect, 2)
            varConnect(i, j) = -1;
        end
    end
    
    % 完成其他信息
    for i = 1:varNum
        cnt = 1;
        while cnt <= size(varConnect, 2) && varConnect(i, cnt) ~= -1
            chk_idx = varConnect(i, cnt) + 1;  % 转换为MATLAB索引
            if chk_idx > 0 && chk_idx <= chkNum
                degChk(chk_idx) = degChk(chk_idx) + 1;
                if degChk(chk_idx) <= size(chkConnect, 2)
                    chkConnect(chk_idx, degChk(chk_idx)) = i - 1;  % 保持0-based值
                    varConnectPosition(i, cnt) = degChk(chk_idx) - 1;
                    chkConnectPosition(chk_idx, degChk(chk_idx)) = cnt - 1;
                end
                degVar(i) = degVar(i) + 1;
            end
            cnt = cnt + 1;
        end
    end
end

% 计算空间耦合正则LDPC的校验节点数
function chkNum = chkSize_SC_regLDPC(l, m, L)
    if mod(l, 2) == 1
        l1 = (l - 1) / 2;
        l2 = (l - 1) / 2;
    else
        l1 = (l - 2) / 2;
        l2 = l / 2;
    end
    chkNum = m * (L + l1 + l2);
end

% AWGN信道仿真
function [received_signal, llr] = AWGN_channel(transmitted_bits, snr_linear)
    n = length(transmitted_bits);
    
    % BPSK调制: 0 -> +1, 1 -> -1
    modulated_signal = 1 - 2 * transmitted_bits;
    
    % 计算噪声功率
    noise_power = 1 / snr_linear;
    
    % 添加AWGN噪声
    noise = sqrt(noise_power) * randn(n, 1);
    received_signal = modulated_signal + noise;
    
    % 计算LLR (对数似然比)
    % LLR = 2 * received_signal / noise_power
    llr = 2 * received_signal / noise_power;
end

% BP解码器 (AWGN版本)
function decodedword = BPDecoder_AWGN(...
    llr, varConnect, chkConnect, varConnectPosition, chkConnectPosition, ...
    degVar, degChk, varNum, chkNum, memorySizeVar, memorySizeChk, MAX_BPIter, MAX_UnchangedIter)
    
    % 初始化消息
    messageVar2Chk = zeros(varNum, memorySizeVar);
    messageChk2Var = zeros(chkNum, memorySizeChk);
    prevword = zeros(varNum, 1);
    
    iter = 0;
    unchangedCnt = 0;
    
    while iter < MAX_BPIter && unchangedCnt < MAX_UnchangedIter
        iter = iter + 1;
        
        % 变量节点到校验节点的消息
        for i = 1:varNum
            for j = 1:degVar(i)
                chk_idx = varConnect(i, j) + 1;
                pos_idx = varConnectPosition(i, j) + 1;
                if chk_idx > 0 && chk_idx <= chkNum && pos_idx > 0 && pos_idx <= memorySizeChk
                    % 计算来自其他校验节点的消息总和
                    sumMessage = llr(i);
                    for k = 1:degVar(i)
                        if k ~= j
                            other_chk_idx = varConnect(i, k) + 1;
                            other_pos_idx = varConnectPosition(i, k) + 1;
                            if other_chk_idx > 0 && other_chk_idx <= chkNum && other_pos_idx > 0 && other_pos_idx <= memorySizeChk
                                sumMessage = sumMessage + messageChk2Var(other_chk_idx, other_pos_idx);
                            end
                        end
                    end
                    messageVar2Chk(i, j) = sumMessage;
                end
            end
        end
        
        % 校验节点到变量节点的消息
        for a = 1:chkNum
            for j = 1:degChk(a)
                var_idx = chkConnect(a, j) + 1;
                pos_idx = chkConnectPosition(a, j) + 1;
                if var_idx > 0 && var_idx <= varNum && pos_idx > 0 && pos_idx <= memorySizeVar
                    % 计算来自其他变量节点的消息
                    sign_product = 1;
                    min_abs = inf;
                    
                    for k = 1:degChk(a)
                        if k ~= j
                            other_var_idx = chkConnect(a, k) + 1;
                            other_pos_idx = chkConnectPosition(a, k) + 1;
                            if other_var_idx > 0 && other_var_idx <= varNum && other_pos_idx > 0 && other_pos_idx <= memorySizeVar
                                msg = messageVar2Chk(other_var_idx, other_pos_idx);
                                sign_product = sign_product * sign(msg);
                                min_abs = min(min_abs, abs(msg));
                            end
                        end
                    end
                    
                    % 使用min-sum近似
                    messageChk2Var(a, j) = sign_product * min_abs;
                end
            end
        end
        
        % 解码
        decodedword = zeros(varNum, 1);
        
        for i = 1:varNum
            % 计算后验LLR
            posterior_llr = llr(i);
            for j = 1:degVar(i)
                chk_idx = varConnect(i, j) + 1;
                pos_idx = varConnectPosition(i, j) + 1;
                if chk_idx > 0 && chk_idx <= chkNum && pos_idx > 0 && pos_idx <= memorySizeChk
                    posterior_llr = posterior_llr + messageChk2Var(chk_idx, pos_idx);
                end
            end
            
            % 硬判决
            if posterior_llr >= 0
                decodedword(i) = 0;
            else
                decodedword(i) = 1;
            end
        end
        
        % 检查收敛
        changed = ~isequal(decodedword, prevword);
        if ~changed
            unchangedCnt = unchangedCnt + 1;
        else
            unchangedCnt = 0;
            prevword = decodedword;
        end
    end
end

% 保存结果到文件 (SNR版本)
function save_results_snr(filename, snr_vals, BER, FER)
    fid = fopen(filename, 'w');
    if fid == -1
        fprintf('无法打开文件 %s\n', filename);
        return;
    end
    
    fprintf(fid, 'SNR (dB):\n');
    for i = 1:length(snr_vals)
        fprintf(fid, '%.1f\n', snr_vals(i));
    end
    
    fprintf(fid, '\nBER:\n');
    for i = 1:length(BER)
        fprintf(fid, '%.10f\n', BER(i));
    end
    
    fprintf(fid, '\nFER:\n');
    for i = 1:length(FER)
        fprintf(fid, '%.10f\n', FER(i));
    end
    
    fclose(fid);
    fprintf('结果已保存到 %s\n', filename);
end

% 绘制结果 (SNR版本)
function plot_results_snr(snr_vals, BER, FER)
    figure('Position', [100, 100, 1200, 400]);
    
    % BER曲线
    subplot(1, 2, 1);
    semilogy(snr_vals, BER, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('信噪比 SNR (dB)');
    ylabel('误比特率 (BER)');
    title('空间耦合LDPC - BER性能 (AWGN信道)');
    grid on;
    
    % FER曲线
    subplot(1, 2, 2);
    semilogy(snr_vals, FER, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('信噪比 SNR (dB)');
    ylabel('帧错误率 (FER)');
    title('空间耦合LDPC - FER性能 (AWGN信道)');
    grid on;
    
    % 保存图形
    savefig('spatial_coupling_ldpc_snr_results.fig');
    print('-dpng', '-r300', 'spatial_coupling_ldpc_snr_results.png');
    fprintf('图形已保存\n');
end
