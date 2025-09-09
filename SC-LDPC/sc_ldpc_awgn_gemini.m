%****************************************************************************************
%* *
%* EE 388 - Spatial Coupling LDPC (MATLAB Version for AWGN Channel)            *
%* *
%****************************************************************************************

function spatial_coupling_ldpc_awgn_simulation()
    % 参数定义
    m_PROTOGRAPH = 100;
    L_SPATIALCOUPLING = 100;
    k_LDPC_reg_ensemble = 12;
    l_LDPC_reg_ensemble = 2;    % k/l = q, 必须是整数
    
    MAX_BP_ITER = 50;
    MAX_BP_UNCHANGED_ITER = 10;
    SIMULATION_NUM = 2; % 为获得平滑曲线，建议增加此值
    
    % AWGN信道仿真参数 - 信噪比 (dB)
    SNR_dB_vals = [-1,-0.5,0,0.5,1,1.5,2,2.5,3,4, 5,6,7,8];
    
    fprintf('开始空间耦合LDPC在AWGN信道下的仿真...\n');
    
    % 代码构造
    m = m_PROTOGRAPH;
    l = l_LDPC_reg_ensemble;
    k = k_LDPC_reg_ensemble;
    L = L_SPATIALCOUPLING;
    q = k / l;
    
    fprintf('参数: m=%d, l=%d, k=%d, L=%d, q=%d\n', m, l, k, L, q);
    
    varNum = q * L * m;
    chkNum = chkSize_SC_regLDPC(l, m, L);
    
    % 计算码率
    code_rate = (varNum - chkNum) / varNum;
    fprintf('变量节点数: %d, 校验节点数: %d, 码率 R = %.4f\n', varNum, chkNum, code_rate);
    
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
    BER = zeros(1, length(SNR_dB_vals));
    
    % 主仿真循环 (使用 parfor 进行并行计算)
    parfor cnt_snr = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(cnt_snr);
        noiseVar = 1 / (10^(SNR_dB / 10)); % 信号功率为1
        
        total_errors = 0;
        
        for cnt_simul = 1:SIMULATION_NUM
            % 信道仿真
            % 假设发送全零码字，BPSK调制后为全+1信号
            tx_signal = ones(varNum, 1);
            
            % AWGN信道：直接添加高斯噪声
             y = awgn(tx_signal, SNR_dB, 'measured'); % 需要Communications Toolbox
            % 手动实现 AWGN，无需工具箱
            % noise = sqrt(noiseVar) * randn(varNum, 1);
            % y = tx_signal + noise;

            % 计算初始LLR
            Lch = 2 * y / noiseVar;
            
            % BP解码
            decoded_bits = BPDecoder_AWGN(...
                Lch, varConnect, chkConnect, varConnectPosition, chkConnectPosition, ...
                degVar, degChk, varNum, chkNum, l+1, k+1, MAX_BP_ITER, MAX_BP_UNCHANGED_ITER);
            
            % 计算BER (发送的是全零码字)
            error_ber = sum(decoded_bits ~= 0);
            total_errors = total_errors + error_ber;
        end
        
        BER(cnt_snr) = total_errors / (varNum * SIMULATION_NUM);
        fprintf('SNR = %.2f dB, BER = %.10f\n', SNR_dB, BER(cnt_snr));
    end
    
    % 显示最终结果
    fprintf('\n最终结果:\n');
    for cnt_snr = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(cnt_snr);
        fprintf('SNR = %.2f dB, BER = %.10f\n', SNR_dB, BER(cnt_snr));
    end
    
    % 保存结果到文件
    filename = sprintf('Result_AWGN_%d.txt', L_SPATIALCOUPLING);
    save_results_awgn(filename, SNR_dB_vals, BER);
    
    % 绘图
    plot_results_awgn(SNR_dB_vals, BER);
    
    fprintf('仿真完成!\n');
end

%% BP解码器 (AWGN信道 - LLR和积算法)
function decodedword = BPDecoder_AWGN(...
    Lch, varConnect, chkConnect, varConnectPosition, chkConnectPosition, ...
    degVar, degChk, varNum, chkNum, memorySizeVar, memorySizeChk, MAX_BPIter, MAX_UnchangedIter)
    
    % 初始化消息
    messageChk2Var = zeros(chkNum, memorySizeChk); % C2V LLRs
    prev_decodedword = 2 * ones(varNum, 1); % 初始化为与0,1都不同的值
    
    iter = 0;
    unchangedCnt = 0;
    
    % 创建一个与 varConnect 大小相同的矩阵来存储 V2C 消息
    messageVar2Chk = zeros(varNum, memorySizeVar);

    while iter < MAX_BPIter && unchangedCnt < MAX_UnchangedIter
        iter = iter + 1;
        
        % --- 步骤 1: 变量节点到校验节点的消息 (V2C) ---
        for i = 1:varNum
            for j = 1:degVar(i)
                % 计算除了当前边之外的所有C2V消息之和
                c2v_sum = 0;
                for k = 1:degVar(i)
                    if k ~= j
                        chk_idx_k = varConnect(i, k) + 1;
                        pos_idx_k = varConnectPosition(i, k) + 1;
                        c2v_sum = c2v_sum + messageChk2Var(chk_idx_k, pos_idx_k);
                    end
                end
                messageVar2Chk(i, j) = Lch(i) + c2v_sum;
            end
        end
        
        % --- 步骤 2: 校验节点到变量节点的消息 (C2V) ---
        for a = 1:chkNum
            for j = 1:degChk(a)
                % 计算除了当前边之外的所有V2C消息的tanh乘积
                tanh_prod = 1.0;
                for k = 1:degChk(a)
                    if k ~= j
                        var_idx_k = chkConnect(a, k) + 1;
                        pos_idx_k = chkConnectPosition(a, k) + 1;
                        v2c_msg = messageVar2Chk(var_idx_k, pos_idx_k);
                        tanh_prod = tanh_prod * tanh(v2c_msg / 2);
                    end
                end
                
                % 防止 atanh(1) 或 atanh(-1) 导致 Inf
                if abs(tanh_prod) > 0.9999999999
                    tanh_prod = sign(tanh_prod) * 0.9999999999;
                end
                messageChk2Var(a, j) = 2 * atanh(tanh_prod);
            end
        end
        
        % --- 步骤 3: 解码与收敛检查 ---
        L_posterior = zeros(varNum, 1);
        for i = 1:varNum
            c2v_sum_total = 0;
            for j = 1:degVar(i)
                chk_idx = varConnect(i, j) + 1;
                pos_idx = varConnectPosition(i, j) + 1;
                c2v_sum_total = c2v_sum_total + messageChk2Var(chk_idx, pos_idx);
            end
            L_posterior(i) = Lch(i) + c2v_sum_total;
        end
        
        % 硬判决
        decodedword = (L_posterior < 0);
        
        % 检查收敛
        if isequal(decodedword, prev_decodedword)
            unchangedCnt = unchangedCnt + 1;
        else
            unchangedCnt = 0;
            prev_decodedword = decodedword;
        end
    end
end


%% 辅助函数 (部分已修改以适应AWGN仿真)

function save_results_awgn(filename, SNR_dB_vals, BER)
    fid = fopen(filename, 'w');
    if fid == -1
        fprintf('无法打开文件 %s\n', filename);
        return;
    end
    
    fprintf(fid, 'SNR_dB:\n');
    for i = 1:length(SNR_dB_vals)
        fprintf(fid, '%.10f\n', SNR_dB_vals(i));
    end
    
    fprintf(fid, '\nBER:\n');
    for i = 1:length(BER)
        fprintf(fid, '%.10f\n', BER(i));
    end
    
    fclose(fid);
    fprintf('结果已保存到 %s\n', filename);
end

function plot_results_awgn(SNR_dB_vals, BER)
    figure('Position', [100, 100, 800, 600]);
    
    semilogy(SNR_dB_vals, BER, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('信噪比 SNR (dB)');
    ylabel('误比特率 (BER)');
    title('空间耦合LDPC - AWGN信道BER性能');
    grid on;
    
    % 保存图形
    savefig('spatial_coupling_ldpc_results_awgn.fig');
    print('-dpng', '-r300', 'spatial_coupling_ldpc_results_awgn.png');
    fprintf('图形已保存\n');
end

% --- 以下是原始程序中无需修改的辅助函数 ---

function vector = randperm_custom(n)
    vector = 0:(n-1);
    for k = 1:n
        j = randi([k, n]);
        temp = vector(j);
        vector(j) = vector(k);
        vector(k) = temp;
    end
end

function matrix = create_connection_matrix(rows, cols)
    matrix = -ones(rows, cols, 'int32');
end

function varConnect = coupling_regular(l, q, L, varConnect)
    if mod(l, 2) == 1
        l1 = (l - 1) / 2;
        l2 = (l - 1) / 2;
    else
        l1 = (l - 2) / 2;
        l2 = l / 2;
    end
    
    for i = 1:L
        cnt = 1;
        for a = i:(i + l - 1)
            for r = 1:q
                row_idx = (r-1)*L + i;
                if row_idx <= size(varConnect, 1) && cnt <= size(varConnect, 2)
                    varConnect(row_idx, cnt) = a - 1;
                end
            end
            cnt = cnt + 1;
        end
    end
end

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

function [varConnect, chkConnect, degVar, degChk, varConnectPosition, chkConnectPosition] = ...
    complete_graph(varNum, chkNum, k, l, varConnect)
    
    chkConnect = create_connection_matrix(chkNum, k+1);
    varConnectPosition = create_connection_matrix(varNum, l+1);
    chkConnectPosition = create_connection_matrix(chkNum, k+1);
    degVar = zeros(varNum, 1);
    degChk = zeros(chkNum, 1);
    
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
        
        for j = 1:length(unique_connections)
            if j <= size(varConnect, 2)
                varConnect(i, j) = unique_connections(j);
            end
        end
        for j = (length(unique_connections)+1):size(varConnect, 2)
            varConnect(i, j) = -1;
        end
    end
    
    for i = 1:varNum
        cnt = 1;
        while cnt <= size(varConnect, 2) && varConnect(i, cnt) ~= -1
            chk_idx = varConnect(i, cnt) + 1;
            if chk_idx > 0 && chk_idx <= chkNum
                degChk(chk_idx) = degChk(chk_idx) + 1;
                if degChk(chk_idx) <= size(chkConnect, 2)
                    chkConnect(chk_idx, degChk(chk_idx)) = i - 1;
                    varConnectPosition(i, cnt) = degChk(chk_idx) - 1;
                    chkConnectPosition(chk_idx, degChk(chk_idx)) = cnt - 1;
                end
                degVar(i) = degVar(i) + 1;
            end
            cnt = cnt + 1;
        end
    end
end

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