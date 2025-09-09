%****************************************************************************************
%* *
%* EE 388 - Spatial Coupling LDPC (MATLAB Version for AWGN Channel)          *
%* *
%****************************************************************************************

function spatial_coupling_ldpc_awgn_simulation()
    % 参数定义
    m_PROTOGRAPH = 100;
    L_SPATIALCOUPLING = 500;
    k_LDPC_reg_ensemble = 10;
    l_LDPC_reg_ensemble = 5;    % k/l = q, 必须是整数
    
    MAX_BP_ITER = 100;
    % 对于AWGN信道，收敛检查可以简化或移除，通常迭代固定次数
    
    SIMULATION_NUM = 20; % 为了快速演示，设置一个较小值。实际仿真建议设为100或更高。
    
    % ======================== MODIFICATION START ========================
    % 修改: 将仿真参数从擦除概率eps改为信噪比Es/N0 (dB)
    EsN0_dB_vals = [0,2,4,6,8,10,12];
    % ========================= MODIFICATION END =========================
    
    fprintf('开始空间耦合LDPC在AWGN信道下的仿真...\n');
    
    % 代码构造 (这部分与信道无关，保持不变)
    m = m_PROTOGRAPH;
    l = l_LDPC_reg_ensemble;
    k = k_LDPC_reg_ensemble;
    L = L_SPATIALCOUPLING;
    q = k / l;
    
    fprintf('参数: m=%d, l=%d, k=%d, L=%d, q=%d\n', m, l, k, L, q);
    
    varNum = q * L * m;
    chkNum = chkSize_SC_regLDPC(l, m, L);
    
    % 码率计算
    code_rate = (varNum - chkNum) / varNum;
    fprintf('变量节点数: %d, 校验节点数: %d, 码率 R: %.4f\n', varNum, chkNum, code_rate);
    
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
    BER = zeros(1, length(EsN0_dB_vals));
    
    % 主仿真循环
    parfor cnt_esn0 = 1:length(EsN0_dB_vals)
        EsN0_dB = EsN0_dB_vals(cnt_esn0);
        EsN0_linear = 10^(EsN0_dB / 10);
        
        % 噪声方差 sigma^2 = N0/2. 假设 BPSK 信号能量 Es=1, 则 N0 = 1/EsN0_linear.
        sigma2 = 1 / (2 * EsN0_linear);
        
        num_errors = 0;
        
        for cnt_simul = 1:SIMULATION_NUM
            % ======================== MODIFICATION START ========================
            % 修改: 信道仿真 - 使用AWGN信道并计算LLR
            % 假设发送全零码字 (bits = 0), BPSK调制后为全+1信号
            tx_signal = ones(varNum, 1); 
            
            % AWGN信道并计算初始LLR
            [y, LLR_ch] = AWGN_channel(tx_signal, sigma2);
            
            % 修改: 使用适用于AWGN的BP解码器 (Sum-Product Algorithm)
            decoded_bits = BPDecoder_AWGN(LLR_ch, varConnect, chkConnect, ...
                degVar, degChk, varNum, chkNum, l+1, k+1, MAX_BP_ITER);
            
            % 计算BER (发送的是全零码字)
            error_ber = sum(decoded_bits ~= 0);
            num_errors = num_errors + error_ber;
            % ========================= MODIFICATION END =========================
         end
        
        BER(cnt_esn0) = num_errors / (varNum * SIMULATION_NUM);
        fprintf('Es/N0 = %.2f dB, BER = %.10f\n', EsN0_dB, BER(cnt_esn0));
    end
    
    % 显示最终结果
    fprintf('\n最终结果:\n');
    for cnt_esn0 = 1:length(EsN0_dB_vals)
        fprintf('Es/N0 = %.2f dB, BER = %.10f\n', EsN0_dB_vals(cnt_esn0), BER(cnt_esn0));
    end
    
    % 保存结果到文件
    filename = sprintf('Result_AWGN_%d.txt', L_SPATIALCOUPLING);
    save_results_awgn(filename, EsN0_dB_vals, BER);
    
    % 绘图
    plot_results_awgn(EsN0_dB_vals, BER);
    
    fprintf('仿真完成!\n');
end

%% 图构造辅助函数 (无需修改)

% 创建Fisher-Yates随机排列
function vector = randperm_custom(n)
    vector = 0:(n-1);
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
    
    chkConnect = create_connection_matrix(chkNum, k+1);
    varConnectPosition = create_connection_matrix(varNum, l+1);
    chkConnectPosition = create_connection_matrix(chkNum, k+1);
    degVar = zeros(varNum, 1, 'int32');
    degChk = zeros(chkNum, 1, 'int32');
    
    for i = 1:varNum
        unique_connections = unique(varConnect(i, varConnect(i,:) ~= -1));
        varConnect(i, :) = -1;
        varConnect(i, 1:length(unique_connections)) = unique_connections;
    end
    
    for i = 1:varNum
        for cnt = 1:size(varConnect, 2)
            if varConnect(i, cnt) == -1, break; end
            
            chk_idx = varConnect(i, cnt) + 1;
            if chk_idx > 0 && chk_idx <= chkNum
                degVar(i) = degVar(i) + 1;
                degChk(chk_idx) = degChk(chk_idx) + 1;
                
                if degChk(chk_idx) <= size(chkConnect, 2)
                    chkConnect(chk_idx, degChk(chk_idx)) = i - 1;
                    varConnectPosition(i, degVar(i)) = degChk(chk_idx) -1;
                    chkConnectPosition(chk_idx, degChk(chk_idx)) = degVar(i) -1;
                end
            end
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

%% ======================== MODIFICATION START ========================
%  新增/修改的函数
% ====================================================================

% 函数: AWGN_channel
% 功能: 模拟BPSK调制、AWGN信道、并计算初始LLR
function [rx_signal, LLR_ch] = AWGN_channel(tx_signal, sigma2)
    % tx_signal: BPSK调制信号 (+1 for bit 0, -1 for bit 1)
    % sigma2: 噪声方差
    noise = sqrt(sigma2) * randn(size(tx_signal));
    rx_signal = tx_signal + noise;
    
    % 计算信道LLR
    % LLR(i) = log(P(x_i=0|y_i)/P(x_i=1|y_i)) = 2*y_i/sigma^2
    LLR_ch = 2 * rx_signal / sigma2;
end

% 函数: BPDecoder_AWGN
% 功能: 适用于AWGN信道的和积(Sum-Product)解码算法
function decoded_bits = BPDecoder_AWGN(LLR_ch, varConnect, chkConnect, ...
    degVar, degChk, varNum, chkNum, memorySizeVar, memorySizeChk, MAX_BPIter)
    
    % 初始化消息
    % messageVar2Chk: 变量节点 -> 校验节点 的消息 (LLR)
    % messageChk2Var: 校验节点 -> 变量节点 的消息 (LLR)
    messageVar2Chk = zeros(varNum, memorySizeVar);
    messageChk2Var = zeros(chkNum, memorySizeChk);
    
    % 初始化: 变量节点广播信道LLR给所有连接的校验节点
    for i = 1:varNum
        for j = 1:degVar(i)
            messageVar2Chk(i, j) = LLR_ch(i);
        end
    end

    % BP迭代
    for iter = 1:MAX_BPIter
        % --- 步骤1: 校验节点更新 ---
        % 从所有连接的变量节点收集消息，并计算传出的消息
        for a = 1:chkNum
            for j = 1:degChk(a)
                var_idx = chkConnect(a, j) + 1; % 目标变量节点
                
                % 计算传给var_idx的消息
                tanh_prod = 1.0;
                for l_neighbor = 1:degChk(a)
                    if l_neighbor == j, continue; end
                    
                    var_neighbor_idx = chkConnect(a, l_neighbor) + 1;
                    
                    % 找到var_neighbor_idx传给校验节点a的消息
                    % 这需要反向查找varConnect
                    pos_in_var = find(varConnect(var_neighbor_idx, 1:degVar(var_neighbor_idx)) == a-1, 1);
                    incoming_msg = messageVar2Chk(var_neighbor_idx, pos_in_var);
                    
                    % 使用tanh规则，并增加数值稳定性处理
                    tanh_val = tanh(incoming_msg / 2);
                    if abs(tanh_val) > 0.99999999
                        tanh_val = sign(tanh_val);
                    end
                    tanh_prod = tanh_prod * tanh_val;
                end
                
                % atanh可能产生Inf，进行裁剪
                if abs(tanh_prod) >= 1.0
                    messageChk2Var(a,j) = sign(tanh_prod) * 30; % 用一个大的LLR值近似
                else
                    messageChk2Var(a,j) = 2 * atanh(tanh_prod);
                end
            end
        end

        % --- 步骤2: 变量节点更新 ---
        % 汇总来自校验节点的消息和信道LLR，计算传出的消息
        for i = 1:varNum
            for j = 1:degVar(i)
                chk_idx = varConnect(i, j) + 1; % 目标校验节点
                
                % 找到校验节点chk_idx在chkConnect中的位置
                pos_in_chk = find(chkConnect(chk_idx, 1:degChk(chk_idx)) == i-1, 1);
                
                % 计算传给chk_idx的消息
                sum_msg = LLR_ch(i);
                for l_neighbor = 1:degVar(i)
                    if l_neighbor == j, continue; end
                    
                    chk_neighbor_idx = varConnect(i, l_neighbor) + 1;
                    pos_neighbor_in_chk = find(chkConnect(chk_neighbor_idx, 1:degChk(chk_neighbor_idx)) == i-1, 1);
                    
                    sum_msg = sum_msg + messageChk2Var(chk_neighbor_idx, pos_neighbor_in_chk);
                end
                messageVar2Chk(i,j) = sum_msg;
            end
        end
    end
    
    % --- 最终解码 ---
    LLR_posterior = zeros(varNum, 1);
    for i = 1:varNum
        sum_msg = LLR_ch(i);
        for j = 1:degVar(i)
            chk_idx = varConnect(i, j) + 1;
            pos_in_chk = find(chkConnect(chk_idx, 1:degChk(chk_idx)) == i-1, 1);
            sum_msg = sum_msg + messageChk2Var(chk_idx, pos_in_chk);
        end
        LLR_posterior(i) = sum_msg;
    end
    
    % 硬判决: LLR > 0 -> bit 0; LLR < 0 -> bit 1
    decoded_bits = (LLR_posterior < 0);
end

% 函数: save_results_awgn
% 功能: 保存AWGN仿真结果到文件
function save_results_awgn(filename, esn0_vals, BER)
    fid = fopen(filename, 'w');
    if fid == -1
        fprintf('无法打开文件 %s\n', filename);
        return;
    end
    
    fprintf(fid, 'Es/N0 (dB):\n');
    fprintf(fid, '%.10f\n', esn0_vals);
    
    fprintf(fid, '\nBER:\n');
    fprintf(fid, '%.10f\n', BER);
    
    fclose(fid);
    fprintf('结果已保存到 %s\n', filename);
end

% 函数: plot_results_awgn
% 功能: 绘制AWGN下的BER曲线
function plot_results_awgn(esn0_vals, BER)
    figure('Name', 'SC-LDPC Performance in AWGN', 'Position', [100, 100, 800, 600]);
    
    semilogy(esn0_vals, BER, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('Es/N0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('空间耦合LDPC在AWGN信道下的BER性能');
    grid on;
    axis([min(esn0_vals) max(esn0_vals) 1e-6 1]); % 调整坐标轴范围
    
    % 保存图形
    savefig('sc_ldpc_awgn_results.fig');
    print('-dpng', '-r300', 'sc_ldpc_awgn_results.png');
    fprintf('图形已保存\n');
end

% ========================= MODIFICATION END =========================