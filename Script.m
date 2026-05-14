%% import
load("C:\Users\Egor\Desktop\MIPT\2nd_semester\Lab_plasma\plasmaero-lab\raw_data_2026-04-16.mat")

%% settings
% с. 364 и с. 398

fs = 32000;

ch_idx  = [1, 3]; % индексы каналов в data.raw
% 1 - CTA
% 2 - акселерометр
% 3 - surface sensor

ch_gain = [1, 10];  % усиления
n_sig = numel(ch_idx);

win      = hamming(4096);
noverlap = 2048;
nfft     = 1024; % количество точек БПФ

sensor_coord=[331.75 160.4];

n_x    = 21;
n_y    = 45;
n_pos  = n_x * n_y;
n_freq = nfft/2 + 1;
meas   = 2;

x_coord = 12; % координаты x и z (выбрали близкие к поверхностному датчику)
z_coord = 24; % используются для расчета когерентности

%% вычисление взаимных спектров

% плоские матрицы (n_freq × n_pos) -> reshape в (n_freq × n_x × n_y)
spec_flat = cell(n_sig, n_sig);
for i = 1:n_sig
    for j = 1:n_sig
        spec_flat{i,j} = zeros(n_freq, n_pos);
    end
end

for xi = 1:n_x
    for yi = 1:n_y
        pos = xi + (yi-1)*n_x;
        S = cell(n_sig, 1);
        for k = 1:n_sig
            sig = data.raw(:, ch_idx(k), xi, yi) * ch_gain(k);
            [S{k}, f] = spectrogram(sig, win, noverlap, nfft, fs);
        end
        for i = 1:n_sig
            for j = 1:n_sig
                spec_flat{i,j}(:, pos) = mean(S{i} .* conj(S{j}), 2);
            end
        end
    end
end

data.f = f;

% data.spec{i,j} -> матрица (n_freq × n_x × n_y)
data.spec = cell(n_sig, n_sig);
for i = 1:n_sig
    for j = 1:n_sig
        data.spec{i,j} = reshape(spec_flat{i,j}, n_freq, n_x, n_y);
    end
end

%% new block
range = 262:268;
% 5-15 (250Hz) ; 100-150 (3KHz) ; 300-350 (8KHz)
% в конце презентации таблицу: волновое число и соответствующее значение бета
%ranges = {8:13, 45:48, 153:157, 195:198, 262:268};

%range_tags = {'250hz', '1.4khz', '4.8khz', '6.1khz', '8.2khz'};

figure(); % вот тут сделать гифку, меняя значение в экспоненте (-2 -> менять на разные с шагом в 10 градусов)
viz = real(squeeze(mean(data.spec{2, 1}(range,:,:)*exp(-2*1i*pi), 1))); % надо в районе 5 кгц и 2 кгц

% contourf(data.x, data.z, viz);

p = pcolor(data.x, data.z, viz);
% p = pcolor(viz);
set(p, "LineStyle", "none", "FaceColor", "interp");
colormap jet;
colorbar;
xlabel('x, мм', 'FontSize', 16);
ylabel('z, мм', 'FontSize', 16);
set(gca, 'FontSize', 16);

% unwrap
%% beta спектр по волновому числу (вдоль z)
% 5-15 (250Hz) ; 100-150 (3KHz) ; 300-350 (8KHz)

% интересные диапазоны: около 250гц 1.4кгц, 2.4кгц, 6.1кгц, 4.8кгц, 8.2кгц
ranges     = {8:13, 45:48, 153:157, 195:198, 262:268};

range_tags = {'250hz', '1.4khz', '4.8khz', '6.1khz', '8.2khz'};

stpz     = data.z(1,2)-data.z(1,1);
beta_vect = [-1/stpz/2 : 1/squeeze((data.z(1,end)-data.z(1,1))) : 1/stpz/2];

for r = 1:numel(ranges)
    range = ranges{r};
    beta_spec = mean(mean((fft(data.spec{2, 1}(range,5:10,:), [],3)), 1), 2);

    figure();
    plot(beta_vect, fftshift(squeeze(abs(beta_spec))));
    title(['f= ',num2str(data.f([min(range),max(range)])'/1e3,'%.2f - %.2f'),' kHz'], 'FontSize', 16)
    xlabel('\beta, мм^{-1}', 'FontSize', 16)
    ylabel('|G_\beta|', 'FontSize', 16)
    set(gca, 'FontSize', 16)
    grid on
    grid minor

    saveas(gcf, ['beta_spec_', range_tags{r}, '.png']);
end

%% корреляционная функция
plot(data.z(5,:), mean(viz(5:8,:), 1), "LineWidth", 1); % показать (корелляционная функция)
set(gca, "FontSize", 16);
xlabel("z, mm");
ylabel("C_{xy}");
title("Correlation function");
grid on;

%% графики автоспектров (точка pos = (5,2)) (точка другая, см. в блоке с переменными)

%z_coord = 10;

pxx1 = abs(data.spec{1,1}(:, x_coord, z_coord));
pxx2 = abs(data.spec{2,2}(:, x_coord, z_coord));

figure();
loglog(f, pxx1, 'DisplayName', 'CTA');
hold on;
loglog(f, pxx2, 'DisplayName', 'surface sensor');
hold off;
xlabel('f, Гц', 'FontSize', 16);
ylabel('G_{xx}', 'FontSize', 16);
title('autospectre', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 16);
set(gca, 'FontSize', 16);
grid on;
xlim([1 fs/2]);

% взаимный спектр (точка pos = (5,2))

Sxy = data.spec{1,2}(:, x_coord, z_coord);

pxy = abs(Sxy);
phi = angle(Sxy);

figure();
subplot(2,1,1);
loglog(f, pxy);
xlabel('f, Гц', 'FontSize', 16);
ylabel('|G_{xy}|', 'FontSize', 16);
title('Амплитуда взаимного спектра', 'FontSize', 16);
set(gca, 'FontSize', 16);
grid on;
xlim([1 fs/2]);

subplot(2,1,2);
plot(f, rad2deg(unwrap(phi)));
xlabel('f, Гц', 'FontSize', 16);
ylabel('Фаза, °', 'FontSize', 16);
title('Фаза взаимного спектра', 'FontSize', 16);
set(gca, 'FontSize', 16);
grid on;
xlim([1 fs/2]);
%ylim([-180 180]);

% разобраться до конца со взаимным спектром, построить для разных точек
% и попытаться найти фазовую скорость по наклону графика

%% спектр когерентности


coh = abs(Sxy) ./ sqrt(pxx1 .* pxx2);

idx_1khz = find(f >= 2700, 1);
coh(idx_1khz:end) = smooth(coh(idx_1khz:end), 15);
coh = smooth(coh, 5);

figure();
semilogx(f, coh, 'LineWidth', 1.5);
xlabel('f, Hz', 'FontSize', 16);
ylabel('\gamma', 'FontSize', 16);
title('coherence spectre (CTA — surface sensor)', 'FontSize', 16);
set(gca, 'FontSize', 16);
hold on;
ranges     = {8:13, 44:49, 150:160, 259:271};
range_tags = {'250 Hz', '1.4 kHz', '4.8 kHz', '8.2 kHz'};

for k = 1:numel(ranges)
    f_lo = f(ranges{k}(1));
    f_hi = f(ranges{k}(end));
    patch([f_lo f_hi f_hi f_lo], [0 0 1 1], [1 0.3 0.3], ...
          'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
    %text(sqrt(f_lo*f_hi), 0.92, range_tags{k}, ...
         %'HorizontalAlignment', 'center', 'FontSize', 12);
end
ylim([0 1]);
xlim([1 fs/2]);
grid on;


%% спектр усиления N(f, β) (не работает, не знаю почему)
% data.spec{2,1} = G_21(f,x,z): кросс-спектр поверхностного датчика (фикс.) × CTA

% FFT по z → амплитуда каждой (f,β) моды на каждой x-позиции
G_hat = fftshift(fft(data.spec{2,1}, [], 3), 3);  % (n_freq × n_x × n_z)
A = abs(G_hat);  % ∝ амплитуда волны
A0 = A(:, 1, :) + eps;

% N(f,β) = max по x усиления ln(A/A0)
N_3d = log(A ./ A0);
N = squeeze(max(N_3d, [], 2));% (n_freq × n_z)

% β-вектор (мм⁻¹)
n_z  = size(data.spec{2,1}, 3);
L_z  = data.z(1,end) - data.z(1,1);
beta = fftshift((-floor(n_z/2) : ceil(n_z/2)-1) / L_z);

figure();
pcolor(beta, data.f / 1e3, max(N, 0));
shading interp;
colormap jet;
colorbar;
xlabel('\beta, мм^{-1}', 'FontSize', 16);
ylabel('f, кГц', 'FontSize', 16);
title('Спектр усиления N = ln(A/A_0)', 'FontSize', 16);
set(gca, 'FontSize', 16);
ylim([0 2]);

%% гифка: sweep по фазе с шагом 10 градусов
% формула: exp(-n*1i*pi/18), шаг 10° => n = 0,1,...,35
%ranges = {8:13, 45:48, 153:157, 195:198, 262:268};

%range_tags = {'250hz', '1.4khz', '4.8khz', '6.1khz', '8.2khz'};
range = 262:268;
gif_filename = 'phase_sweep_8.2khz.gif';

% предвычислить все кадры и найти глобальные пределы цветовой шкалы
n_frames = 36;
all_viz = zeros(size(data.spec{2,1}, 2), size(data.spec{2,1}, 3), n_frames);
for n = 0:n_frames-1
    all_viz(:,:,n+1) = real(squeeze(mean(data.spec{2,1}(range,:,:) * exp(-n*1i*pi/18), 1)));
end
clim_val = max(abs(all_viz(:)));

fig = figure();
for n = 0:n_frames-1
    clf;
    p = pcolor(data.x, data.z, all_viz(:,:,n+1));
    set(p, 'LineStyle', 'none', 'FaceColor', 'interp');
    colormap(jet);
    caxis([-clim_val, clim_val]);
    colorbar;
    title(sprintf('\\phi = %d°  (n = %d)', n*10, n), 'FontSize', 16);
    xlabel('x, мм', 'FontSize', 16);
    ylabel('z, мм', 'FontSize', 16);
    set(gca, 'FontSize', 16);
    drawnow;
    drawpoint("Position", sensor_coord, "Color", "w");

    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if n == 0
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end
