%% LDPC与空间耦合LDPC性能对比仿真
% 统一设置信噪比，在同一图中比较两种LDPC的性能
% 支持可调码率、字体大小和线条粗细

function ldpc_comparison_simulation()
    %% 清理环境
    clear all;
    clc;
    
    %% 可调参数设置
    % 信噪比范围 (dB)
    SNR_dB_values = [ 0, 1, 2, 3, 4, 5, 6, 7, 8,9,10,11,12,13,14,15];
    
    % 普通LDPC参数
    LDPC_BLK_SIZE = 256;           % 块大小
    LDPC_CODE_RATE = '5/6';        % 码率 (可调)
    LDPC_NUM_ITER = 100;           % 仿真迭代次数
    
    % 空间耦合LDPC参数
    SC_M_PROTOGRAPH = 100;         % 原型图参数
    SC_L_SPATIALCOUPLING = 100;    % 空间耦合长度
    SC_K_LDPC_REG = 12;            % 变量节点度数
    SC_L_LDPC_REG = 2;             % 校验节点度数 (k/l = q)
    SC_SIMULATION_NUM = 10;         % 仿真次数
    
    % 绘图参数 (可调)
    FONT_SIZE = 14;                % 字体大小
    LINE_WIDTH = 2;                % 线条粗细
    MARKER_SIZE = 8;               % 标记大小
    
    % 添加路径
    if ~isdeployed
        addpath('./codes');
        addpath('./SC-LDPC');
    end
    
    fprintf('开始LDPC与空间耦合LDPC性能对比仿真...\n');
    fprintf('信噪比范围: %.1f 到 %.1f dB\n', min(SNR_dB_values), max(SNR_dB_values));
    fprintf('普通LDPC码率: %s\n', LDPC_CODE_RATE);
    
    %% 运行普通LDPC仿真
    fprintf('\n=== 开始普通LDPC仿真 ===\n');
    [ldpc_snr_values, ldpc_ber_results] = run_regular_ldpc_simulation(...
        LDPC_BLK_SIZE, LDPC_CODE_RATE, LDPC_NUM_ITER, SNR_dB_values);
    
    %% 运行空间耦合LDPC仿真
    fprintf('\n=== 开始空间耦合LDPC仿真 ===\n');
    [sc_snr_values, sc_ber_results] = run_spatial_coupling_ldpc_simulation(...
        SC_M_PROTOGRAPH, SC_L_SPATIALCOUPLING, SC_K_LDPC_REG, SC_L_LDPC_REG, ...
        SC_SIMULATION_NUM, SNR_dB_values);
    
    %% 统一绘图
    fprintf('\n=== 绘制对比图 ===\n');
    plot_comparison_results(ldpc_snr_values, ldpc_ber_results, ...
                           sc_snr_values, sc_ber_results, ...
                           LDPC_CODE_RATE, FONT_SIZE, LINE_WIDTH, MARKER_SIZE);
    
    %% 保存结果
    save_comparison_results(ldpc_snr_values, ldpc_ber_results, ...
                           sc_snr_values, sc_ber_results, LDPC_CODE_RATE);
    
    fprintf('\n仿真完成！\n');
end

%% 普通LDPC仿真函数
function [snr_values, ber_results] = run_regular_ldpc_simulation(blkSize, codeRate, numIter, snr_values)
    % 星座大小
    M = 4;
    
    % 获取LDPC结构
    LDPC = ldpcGet(blkSize, codeRate);
    
    % 初始化结果数组
    ber_results = zeros(size(snr_values));
    
    % 对每个信噪比值进行仿真
    parfor snr_idx = 1:length(snr_values)
        snr_dB = snr_values(snr_idx);
        
        numErr = 0;
        
        fprintf('普通LDPC: 仿真 SNR = %.1f dB...\n', snr_dB);
        
        % 仿真循环
        for i = 1:numIter
            % 生成随机数据
            data = randi([0 1], 1, LDPC.numInfBits);
            
            % 编码
            dataEnc = ldpcEncode(data, LDPC);
            
            % QAM映射
            dataMod = qammod(dataEnc(:), M, 'InputType', 'bit', 'UnitAveragePower', true);
            
            % AWGN
            dataRx = awgn(dataMod, snr_dB);
            
            % LLR解映射
            dataLlr = qamdemod(dataRx, M, 'OutputType', 'llr', 'UnitAveragePower', true);
            
            % 解码
            dataHat = ldpcDecode(dataLlr', LDPC);
            
            % 计算误码数
            numErr = numErr + sum(abs(dataHat - data));
        end
        
        % 计算BER
        ber = numErr / (numIter * LDPC.numInfBits);
        
        % 如果误码率为0，设定为1e-7
        if ber == 0
            ber = 1e-7;
        end
        
        ber_results(snr_idx) = ber;
        fprintf('普通LDPC: SNR = %.1f dB, BER = %.2e\n', snr_dB, ber);
    end
end

%% 空间耦合LDPC仿真函数
function [snr_values, ber_results] = run_spatial_coupling_ldpc_simulation(...
    m_protograph, L_spatialcoupling, k_ldpc_reg, l_ldpc_reg, ...
    simulation_num, snr_dB_vals)
    
    % 参数计算
    m = m_protograph;
    l = l_ldpc_reg;
    k = k_ldpc_reg;
    L = L_spatialcoupling;
    q = k / l;
    
    fprintf('空间耦合LDPC参数: m=%d, l=%d, k=%d, L=%d, q=%d\n', m, l, k, L, q);
    
    varNum = q * L * m;
    chkNum = chkSize_SC_regLDPC(l, m, L);
    
    % 计算码率
    code_rate = (varNum - chkNum) / varNum;
    fprintf('空间耦合LDPC码率 R = %.4f\n', code_rate);
    
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
    
    fprintf('空间耦合LDPC代码构造完成.\n');
    
    % 仿真参数
    MAX_BP_ITER = 50;
    MAX_BP_UNCHANGED_ITER = 10;
    
    % 初始化结果
    snr_values = snr_dB_vals;
    ber_results = zeros(1, length(snr_dB_vals));
    
    % 主仿真循环
    parfor cnt_snr = 1:length(snr_dB_vals)
        SNR_dB = snr_dB_vals(cnt_snr);
        noiseVar = 1 / (10^(SNR_dB / 10));
        
        total_errors = 0;
        
        fprintf('空间耦合LDPC: 仿真 SNR = %.1f dB...\n', SNR_dB);
        
        for cnt_simul = 1:simulation_num
            % 发送全零码字，BPSK调制后为全+1信号
            tx_signal = ones(varNum, 1);
            
            % AWGN信道
            y = awgn(tx_signal, SNR_dB, 'measured');
            
            % 计算初始LLR
            Lch = 2 * y / noiseVar;
            
            % BP解码
            decoded_bits = BPDecoder_AWGN(...
                Lch, varConnect, chkConnect, varConnectPosition, chkConnectPosition, ...
                degVar, degChk, varNum, chkNum, l+1, k+1, MAX_BP_ITER, MAX_BP_UNCHANGED_ITER);
            
            % 计算BER
            error_ber = sum(decoded_bits ~= 0);
            total_errors = total_errors + error_ber;
        end
        
        ber_results(cnt_snr) = total_errors / (varNum * simulation_num);
        
        % 如果误码率为0，设定为1e-7
        if ber_results(cnt_snr) == 0
            ber_results(cnt_snr) = 1e-7;
        end
        
        fprintf('空间耦合LDPC: SNR = %.1f dB, BER = %.2e\n', SNR_dB, ber_results(cnt_snr));
    end
end

%% 对比绘图函数
function plot_comparison_results(ldpc_snr, ldpc_ber, sc_snr, sc_ber, code_rate, font_size, line_width, marker_size)
    % 创建图形窗口
    figure('Position', [100, 100, 1000, 700]);
    
    % 处理误码率为0的情况，将其设置为1e-7用于绘图
    ldpc_ber_plot = ldpc_ber;
    ldpc_ber_plot(ldpc_ber == 0) = 1e-7;
    
    sc_ber_plot = sc_ber;
    sc_ber_plot(sc_ber == 0) = 1e-7;
    
    % 绘制普通LDPC结果
    semilogy(ldpc_snr, ldpc_ber_plot, 'bo-', 'LineWidth', line_width, 'MarkerSize', marker_size, ...
             'DisplayName', sprintf('普通LDPC (码率 %s)', code_rate));
    hold on;
    
    % 绘制空间耦合LDPC结果
    semilogy(sc_snr, sc_ber_plot, 'rs-', 'LineWidth', line_width, 'MarkerSize', marker_size, ...
             'DisplayName', sprintf('空间耦合LDPC(码率 %s)',code_rate));
    
    % 设置图形属性
    grid on;
    xlabel('信噪比 (dB)', 'FontSize', font_size, 'FontWeight', 'bold');
    ylabel('误比特率 (BER)', 'FontSize', font_size, 'FontWeight', 'bold');
    title(sprintf('LDPC与空间耦合LDPC性能对比 (码率 %s)', code_rate), ...
          'FontSize', font_size+2, 'FontWeight', 'bold');
    
    % 设置图例
    legend('Location', 'best', 'FontSize', font_size, 'FontWeight', 'bold');
    
    % 设置坐标轴
    xlim([min([ldpc_snr, sc_snr])-0.5, max([ldpc_snr, sc_snr])+0.5]);
    ylim([1e-8, 1]);
    
    % 设置坐标轴字体
    set(gca, 'FontSize', font_size-2, 'FontWeight', 'bold');
    
    % 保存图形
    savefig('ldpc_comparison_results.fig');
    print('-dpng', '-r300', 'ldpc_comparison_results.png');
    fprintf('对比图已保存为 ldpc_comparison_results.fig 和 ldpc_comparison_results.png\n');
end

%% 保存结果函数
function save_comparison_results(ldpc_snr, ldpc_ber, sc_snr, sc_ber, code_rate)
    % 保存到文本文件
    filename = 'ldpc_comparison_results.txt';
    fid = fopen(filename, 'w');
    
    if fid == -1
        fprintf('无法打开文件 %s\n', filename);
        return;
    end
    
    fprintf(fid, 'LDPC与空间耦合LDPC性能对比结果 (码率 %s)\n', code_rate);
    fprintf(fid, '==========================================\n\n');
    
    fprintf(fid, '普通LDPC结果:\n');
    fprintf(fid, 'SNR (dB)\tBER\n');
    for i = 1:length(ldpc_snr)
        fprintf(fid, '%.2f\t\t%.2e\n', ldpc_snr(i), ldpc_ber(i));
    end
    
    fprintf(fid, '\n空间耦合LDPC结果:\n');
    fprintf(fid, 'SNR (dB)\tBER\n');
    for i = 1:length(sc_snr)
        fprintf(fid, '%.2f\t\t%.2e\n', sc_snr(i), sc_ber(i));
    end
    
    fclose(fid);
    fprintf('结果已保存到 %s\n', filename);
end

%% 空间耦合LDPC辅助函数 (从原文件复制)

function decodedword = BPDecoder_AWGN(...
    Lch, varConnect, chkConnect, varConnectPosition, chkConnectPosition, ...
    degVar, degChk, varNum, chkNum, memorySizeVar, memorySizeChk, MAX_BPIter, MAX_UnchangedIter)
    
    % 初始化消息
    messageChk2Var = zeros(chkNum, memorySizeChk);
    prev_decodedword = 2 * ones(varNum, 1);
    
    iter = 0;
    unchangedCnt = 0;
    
    messageVar2Chk = zeros(varNum, memorySizeVar);

    while iter < MAX_BPIter && unchangedCnt < MAX_UnchangedIter
        iter = iter + 1;
        
        % 变量节点到校验节点的消息
        for i = 1:varNum
            for j = 1:degVar(i)
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
        
        % 校验节点到变量节点的消息
        for a = 1:chkNum
            for j = 1:degChk(a)
                tanh_prod = 1.0;
                for k = 1:degChk(a)
                    if k ~= j
                        var_idx_k = chkConnect(a, k) + 1;
                        pos_idx_k = chkConnectPosition(a, k) + 1;
                        v2c_msg = messageVar2Chk(var_idx_k, pos_idx_k);
                        tanh_prod = tanh_prod * tanh(v2c_msg / 2);
                    end
                end
                
                if abs(tanh_prod) > 0.9999999999
                    tanh_prod = sign(tanh_prod) * 0.9999999999;
                end
                messageChk2Var(a, j) = 2 * atanh(tanh_prod);
            end
        end
        
        % 解码与收敛检查
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
        
        decodedword = (L_posterior < 0);
        
        if isequal(decodedword, prev_decodedword)
            unchangedCnt = unchangedCnt + 1;
        else
            unchangedCnt = 0;
            prev_decodedword = decodedword;
        end
    end
end

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
