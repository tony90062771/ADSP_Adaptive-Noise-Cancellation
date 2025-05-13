% 讀取音頻文件
[signal_with_noise, Fs] = audioread('signal_with_noise.wav');
[reference_noise, ~] = audioread('reference_noise.wav');
       
% 初始化變數
filter_order = 10; % 濾波器階數
w_nlms = zeros(filter_order, 1); % NLMS權重向量
w_rls = zeros(filter_order, 1); % RLS權重向量
mu = 0.008; % NLMS步長因子
eps = 0.01; % 計算 NLMS中的步長，避免分母為零的情況
P = eye(filter_order); % RLS的反矩陣
lambda = 0.99; % RLS的遺忘因子

% 儲存權重向量的歷史
w_nlms_history = zeros(filter_order, length(signal_with_noise));
w_rls_history = zeros(filter_order, length(signal_with_noise));

% 適應性濾波
for n = filter_order:length(signal_with_noise)
    % 提取當前的輸入向量
    u = reference_noise(n:-1:n-filter_order+1);
    
    % NLMS
    y = w_nlms' * u;
    e_nlms = signal_with_noise(n) - y;
    w_nlms = w_nlms + (mu/(u'*u + eps)) * e_nlms * u;
    w_nlms_history(:, n) = w_nlms; % 儲存權重向量
    
    % RLS
    k = (P * u) / (lambda + u' * P * u);
    y = w_rls' * u;
    e_rls = signal_with_noise(n) - y;
    w_rls = w_rls + k * e_rls;
    P = (P - k * u' * P) ;                                                                          %/lambda 這樣分析圖才OK
    w_rls_history(:, n) = w_rls; % 儲存權重向量
end

% 繪製權重向量的歷史
figure;
hold on;
colors = jet(filter_order);
for k = 1:filter_order
    plot(w_nlms_history(k, :), 'Color', colors(k, :));
end
title('NLMS Weight History');
xlabel('n');
ylabel('w_n(k)');
legend('w_1', 'w_2', 'w_3', 'w_4', 'w_5', 'w_6', 'w_7', 'w_8', 'w_9', 'w_10');
hold off;

figure;
hold on;
for k = 1:filter_order
    plot(w_rls_history(k, :), 'Color', colors(k, :));
end
title('RLS Weight History');
xlabel('n');
ylabel('w_n(k)');
legend('w_1', 'w_2', 'w_3', 'w_4', 'w_5', 'w_6', 'w_7', 'w_8', 'w_9', 'w_10');
hold off;

% 寫入降噪後的音頻文件
audiowrite('1226_nlms.wav', e_nlms_array, Fs);
audiowrite('1226_rls.wav', e_rls_array, Fs);

% 繪製原始音訊
figure;
subplot(3,1,1);
plot(signal_with_noise);
title('Original Signal');
xlabel('Time (samples)');
ylabel('Amplitude');

% 繪製誤差信號
subplot(3,1,2);
plot(e_nlms_array);
title('NLMS Error Signal');
xlabel('Time (samples)');
ylabel('Amplitude');

subplot(3,1,3);
plot(e_rls_array);
title('RLS Error Signal');
xlabel('Time (samples)');
ylabel('Amplitude');
