%****************************************************************************************
%*                                                                                      *
%*                      EE 388 - Spatial Coupling LDPC (MATLAB Version)               *
%*                                                                                      *
%****************************************************************************************

function spatial_coupling_ldpc_simulation()
    % 参数定义
    m_PROTOGRAPH = 100;
    L_SPATIALCOUPLING = 500;
    k_LDPC_reg_ensemble = 10;
    l_LDPC_reg_ensemble = 5;    % k/l = q, 必须是整数
    
    MAX_BP_ITER = 200;
    MAX_BP_UNCHANGED_ITER = 10;
    SIMULATION_NUM = 20;
    
    eps_vals = [0, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, ...
                0.38, 0.39, 0.4,  0.42, 0.43,  0.45, 0.5, ...
                0.6, 0.8, 1];
    
    fprintf('开始空间耦合LDPC仿真...\n');
    
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
    BER = zeros(1, length(eps_vals));
    EXIT = zeros(1, length(eps_vals));
    
    progress = 0;
    total_sims = SIMULATION_NUM * length(eps_vals);
    
    % 主仿真循环
    parfor cnt_eps = 1:length(eps_vals)
        eps = eps_vals(cnt_eps);
%         fprintf('\n测试 eps = %.3f\n', eps);
        
        for cnt_simul = 1:SIMULATION_NUM
            % 信道仿真 - 假设发送全零码字
            y = BEC_channel(eps, varNum);
            
            % BP解码
            [decodedword, decodedword_without_yi] = BPDecoder_BEC(...
                y, varConnect, chkConnect, varConnectPosition, chkConnectPosition, ...
                degVar, degChk, varNum, chkNum, l+1, k+1, MAX_BP_ITER, MAX_BP_UNCHANGED_ITER);
            
            % 计算BER
            error_ber = sum(decodedword ~= 0);
            BER(cnt_eps) = BER(cnt_eps) + error_ber / varNum;
            
            % 计算EXIT
            error_exit = sum(decodedword_without_yi ~= 0);
            EXIT(cnt_eps) = EXIT(cnt_eps) + error_exit / varNum;
            
%             % 显示进度
%             progress = progress + 1;
%             if mod(progress, 50) == 0  % 每50次仿真显示一次进度
%                 fprintf('完成 %.1f%% (%d/%d)\n', progress*100/total_sims, progress, total_sims);
%             end
         end
%         
        BER(cnt_eps) = BER(cnt_eps) / SIMULATION_NUM;
        EXIT(cnt_eps) = EXIT(cnt_eps) / SIMULATION_NUM;
        fprintf('Eps = %.5f, BER = %.10f, EXIT = %.10f\n', eps, BER(cnt_eps), EXIT(cnt_eps));
    end
    
    % 显示最终结果
    fprintf('\n最终结果:\n');
    for cnt_eps = 1:length(eps_vals)
        eps = eps_vals(cnt_eps);
        fprintf('Eps = %.5f, BER = %.10f, EXIT = %.10f\n', eps, BER(cnt_eps), EXIT(cnt_eps));
    end
    
    % 保存结果到文件
    filename = sprintf('Result_%d.txt', L_SPATIALCOUPLING);
    save_results(filename, eps_vals, BER, EXIT);
    
    % 绘图
    plot_results(eps_vals, BER, EXIT);
    
    fprintf('仿真完成!\n');
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

% BEC信道仿真
function output = BEC_channel(eps, n)
    output = zeros(n, 1);
    erasure_prob = rand(n, 1);
    output(erasure_prob >= (1 - eps)) = -1;  % -1表示擦除
end

% BP解码器
function [decodedword, decodedword_without_yi] = BPDecoder_BEC(...
    y, varConnect, chkConnect, varConnectPosition, chkConnectPosition, ...
    degVar, degChk, varNum, chkNum, memorySizeVar, memorySizeChk, MAX_BPIter, MAX_UnchangedIter)
    
    % 初始化消息
    messageVar2Chk = zeros(varNum, memorySizeVar);
    messageChk2Var = -ones(chkNum, memorySizeChk);
    prevword = 2 * ones(varNum, 1);
    
    iter = 0;
    unchangedCnt = 0;
    
    while iter < MAX_BPIter && unchangedCnt < MAX_UnchangedIter
        iter = iter + 1;
        
        % 变量节点到校验节点的消息
        for i = 1:varNum
            sumMessage = 0;
            for j = 1:degVar(i)
                chk_idx = varConnect(i, j) + 1;
                pos_idx = varConnectPosition(i, j) + 1;
                if chk_idx > 0 && chk_idx <= chkNum && pos_idx > 0 && pos_idx <= memorySizeChk
                    sumMessage = sumMessage + messageChk2Var(chk_idx, pos_idx);
                end
            end
            
            for j = 1:degVar(i)
                chk_idx = varConnect(i, j) + 1;
                pos_idx = varConnectPosition(i, j) + 1;
                if chk_idx > 0 && chk_idx <= chkNum && pos_idx > 0 && pos_idx <= memorySizeChk
                    if (y(i) == 0) || ((sumMessage - messageChk2Var(chk_idx, pos_idx)) ~= (1 - degVar(i)))
                        messageVar2Chk(i, j) = 0;
                    else
                        messageVar2Chk(i, j) = -1;
                    end
                end
            end
        end
        
        % 校验节点到变量节点的消息
        for a = 1:chkNum
            sumMessage = 0;
            for j = 1:degChk(a)
                var_idx = chkConnect(a, j) + 1;
                pos_idx = chkConnectPosition(a, j) + 1;
                if var_idx > 0 && var_idx <= varNum && pos_idx > 0 && pos_idx <= memorySizeVar
                    sumMessage = sumMessage + messageVar2Chk(var_idx, pos_idx);
                end
            end
            
            for j = 1:degChk(a)
                var_idx = chkConnect(a, j) + 1;
                pos_idx = chkConnectPosition(a, j) + 1;
                if var_idx > 0 && var_idx <= varNum && pos_idx > 0 && pos_idx <= memorySizeVar
                    if (sumMessage - messageVar2Chk(var_idx, pos_idx)) ~= 0
                        messageChk2Var(a, j) = -1;
                    else
                        messageChk2Var(a, j) = 0;
                    end
                end
            end
        end
        
        % 解码
        decodedword = zeros(varNum, 1);
        decodedword_without_yi = zeros(varNum, 1);
        
        for i = 1:varNum
            sumMessage = 0;
            for j = 1:degVar(i)
                chk_idx = varConnect(i, j) + 1;
                pos_idx = varConnectPosition(i, j) + 1;
                if chk_idx > 0 && chk_idx <= chkNum && pos_idx > 0 && pos_idx <= memorySizeChk
                    sumMessage = sumMessage + messageChk2Var(chk_idx, pos_idx);
                end
            end
            
            if y(i) == 0 || sumMessage ~= (-degVar(i))
                decodedword(i) = 0;
            else
                decodedword(i) = -1;
            end
            
            if sumMessage ~= (-degVar(i))
                decodedword_without_yi(i) = 0;
            else
                decodedword_without_yi(i) = -1;
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

% 保存结果到文件
function save_results(filename, eps_vals, BER, EXIT)
    fid = fopen(filename, 'w');
    if fid == -1
        fprintf('无法打开文件 %s\n', filename);
        return;
    end
    
    fprintf(fid, 'Eps:\n');
    for i = 1:length(eps_vals)
        fprintf(fid, '%.10f\n', eps_vals(i));
    end
    
    fprintf(fid, '\nEXIT:\n');
    for i = 1:length(EXIT)
        fprintf(fid, '%.10f\n', EXIT(i));
    end
    
    fprintf(fid, '\nBER:\n');
    for i = 1:length(BER)
        fprintf(fid, '%.10f\n', BER(i));
    end
    
    fclose(fid);
    fprintf('结果已保存到 %s\n', filename);
end

% 绘制结果
function plot_results(eps_vals, BER, EXIT)
    figure('Position', [100, 100, 1200, 400]);
    
    % BER曲线
    subplot(1, 2, 1);
    semilogy(eps_vals, BER, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('擦除概率 ε');
    ylabel('误比特率 (BER)');
    title('空间耦合LDPC - BER性能');
    grid on;
    
    % EXIT曲线
    subplot(1, 2, 2);
    plot(eps_vals, EXIT, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('擦除概率 ε');
    ylabel('EXIT函数');
    title('空间耦合LDPC - EXIT函数');
    grid on;
    
    % 保存图形
    savefig('spatial_coupling_ldpc_results.fig');
    print('-dpng', '-r300', 'spatial_coupling_ldpc_results.png');
    fprintf('图形已保存\n');
end