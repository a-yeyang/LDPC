%% Clean up
clear all;
clc;

%% Init
if ~isdeployed
    addpath('./codes');
end

% Constellation size
M = 4;

% LDPC config
blkSize = 256;
codeRate = '5/6';

% Get LDPC struct
LDPC = ldpcGet(blkSize, codeRate);

% Simulation parameters
% ebno_values = [0, 2, 4, 6,8,10,12,14];  % 信噪比数值
ebno_values=[-1,-0.5,0,0.5,1,1.5,2,2.5,3,4, 5,6,7,8];
numIter = 1e2;  % 仿真迭代次数
ber_results = zeros(size(ebno_values));  % 存储误码率结果

%% 对每个信噪比值进行仿真
parfor ebno_idx = 1:length(ebno_values)
    ebno = ebno_values(ebno_idx);
    numErr = 0;
    
    % Convert E_b/N_0 to some SNR
    snr = ebno + 10*log10(log2(M)) + 10*log10(str2num(codeRate));
    
    fprintf('正在仿真 Eb/N0 = %d dB...\n', ebno);
    
    % Simulate
    for i = 1:numIter
        
        % Generate random data
        data = randi([0 1], 1, LDPC.numInfBits);

        % Encode
        dataEnc = ldpcEncode(data, LDPC);

        % QAM mapping
        dataMod = qammod(dataEnc(:), M, 'InputType', 'bit', 'UnitAveragePower', true);

        % AWGN
        dataRx = awgn(dataMod, snr);

        % LLR demapping
        dataLlr = qamdemod(dataRx, M, 'OutputType', 'llr', 'UnitAveragePower', true);

        % Decode
        dataHat = ldpcDecode(dataLlr', LDPC);

        % Count number of bit errors
        numErr = numErr + sum(abs(dataHat - data));
        
    end
    
    % Calculate BER
    ber = numErr / (numIter * LDPC.numInfBits);
    
    % 如果误码率为0，设定为1e-7
    if ber == 0
        ber = 1e-7;
    end
    
    ber_results(ebno_idx) = ber;
    fprintf('Eb/N0 = %d dB, BER = %.2e\n', ebno, ber);
end

%% 绘制误码率曲线
figure;
semilogy(ebno_values, ber_results, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Eb/N0 (dB)');
ylabel('误码率 (BER)');
title('LDPC编码误码率性能曲线');
legend('LDPC (5/6)', 'Location', 'best');

% 设置坐标轴范围
xlim([min(ebno_values)-1, max(ebno_values)+1]);
ylim([1e-8, 1]);

% 显示结果
fprintf('\n仿真结果:\n');
for i = 1:length(ebno_values)
    fprintf('Eb/N0 = %d dB: BER = %.2e\n', ebno_values(i), ber_results(i));
end