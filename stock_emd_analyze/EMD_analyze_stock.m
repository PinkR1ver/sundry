clear all;clc;close all;

% 假设你已经有了一个信号，这里我们创建一个示例信号
% t = 0:0.01:4; % 时间向量
% signal = cos(2*pi*t) + 0.5*cos(4*pi*t) + sin(8*pi*t); % 示例信号

data = readtable('stock.csv');
t = data.Date;
signal = data.Close;

t = 1:length(t);
t = t';

for i = 1:length(signal)
    if isnan(signal(i))
        signal(i) = signal(i - 1);
    end
end


loop_time = 7;

ax = figure("Position", [400 400 800 600]);
subplot(loop_time + 1, 1, 1);
plot(t, signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Original Signal');
grid on;

figure;
plot(t, signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Original Signal');
grid on;

loop_num = 1;

while true


    % get all max value and min to get envelope e_max and e_min

    % get all max value and min value
    maxValue = [];
    maxIndex = [];
    minValue = [];
    minIndex = [];

    for i = 2:length(signal) - 1
        if signal(i) > signal(i - 1) && signal(i) > signal(i + 1)
            maxValue = [maxValue, signal(i)];
            maxIndex = [maxIndex, t(i)];
        end
        
        if signal(i) < signal(i - 1) && signal(i) < signal(i + 1)
            minValue = [minValue, signal(i)];
            minIndex = [minIndex, t(i)];
        end
    end

    % plot the signal and all max value point and min value point

    figure("Position", [400 400 800 600]);
    plot(t, signal, 'LineWidth', 1.5, 'DisplayName', 'Signal');
    hold on;
    plot(maxIndex, maxValue, 'r.', 'MarkerSize', 20, 'DisplayName', 'Max Value');
    plot(minIndex, minValue, 'b.', 'MarkerSize', 20, 'DisplayName', 'Min Value');
    hold off;
    xlabel('Time');
    ylabel('Amplitude');
    title('Signal with Max and Min Value');
    legend('show');
    grid on;

    eMax = interp1(maxIndex, maxValue, t, 'spline');
    eMin = interp1(minIndex, minValue, t, 'spline');

    % plot the signal and envelope e_max and e_min

    figure("Position", [400 400 800 600]);
    plot(t, signal, 'LineWidth', 1.5, 'DisplayName', 'Signal');
    hold on;
    plot(t, eMax, 'r', 'LineWidth', 1.5, 'DisplayName', 'e_{max}');
    plot(t, eMin, 'b', 'LineWidth', 1.5, 'DisplayName', 'e_{min}');
    hold off;
    xlabel('Time');
    ylabel('Amplitude');
    title('Signal with Envelope e_{max} and e_{min}');
    legend('show');
    grid on;

    % calculate the mean of envelope
    m = (eMax + eMin) / 2;

    % plot them

    figure("Position", [400 400 800 600]);
    plot(t, signal, 'LineWidth', 1.5, 'DisplayName', 'Signal');
    hold on;
    plot(t, eMax, 'r', 'LineWidth', 1.5, 'DisplayName', 'e_{max}');
    plot(t, eMin, 'b', 'LineWidth', 1.5, 'DisplayName', 'e_{min}');
    plot(t, m, 'g--', 'LineWidth', 1.5, 'DisplayName', 'm');
    hold off;
    xlabel('Time');
    ylabel('Amplitude');
    title('Signal with Mean m');
    legend('show');
    grid on;

    figure("Position", [400 400 800 600]);
    plot(t, signal, '-.', 'LineWidth', 1.5, 'DisplayName', 'Signal');
    hold on;
    plot(t, m, '--', 'LineWidth', 1.5, 'DisplayName', 'm');
    plot(t, signal - m, 'r-', 'LineWidth', 3, 'DisplayName', 'c(t) = x(t) - m(t)');
    xlabel('Time');
    ylabel('Amplitude');
    title('Signal - m');
    legend('show');
    grid on;

    imf = signal - m;

    figure(ax);
    subplot(loop_time + 1, 1, loop_num + 1);
    plot(t, imf, 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Amplitude');
    title(['IMF', num2str(loop_num)]);
    grid on;

    residuals = signal - imf;
    signal = residuals;

    figure("Position", [400 400 800 600]);
    plot(t, signal, 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Amplitude');
    title(['Residuals', num2str(loop_num)]);
    grid on;

    loop_num = loop_num + 1;

    if loop_num == loop_time
        break;
    end


end

figure(ax);
subplot(loop_time + 1, 1, loop_time + 1);
plot(t, residuals, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Residuals');
grid on;










