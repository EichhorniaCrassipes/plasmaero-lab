%% import
load("C:\Users\Egor\Desktop\MIPT\2nd_semester\Lab_plasma\plasmaero-lab\raw_data_2026-04-16.mat")

%% settings
% с. 364 и с. 398

fs = 32000;

ch_idx  = [1, 3]; % индексы каналов в data.raw
ch_gain = [1, 10];  % усиления
n_sig = numel(ch_idx);

win      = hamming(4096);
noverlap = 2048;
nfft     = 1024; % количество точек БПФ


n_x    = 21;
n_y    = 45;
n_pos  = n_x * n_y;
n_freq = nfft/2 + 1;
meas   = 2;

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
range = 300:350;
% 5-15 (250Hz) ; 100-150 (3KHz) ; 300-350 (8KHz)
% в конце презентации таблицу: волновое число и соответствующее значение бета

figure(); % вот тут сделать гифку, меняя значение в экспоненте (-2 -> менять на разные с шагом в 10 градусов)
viz = real(squeeze(mean(data.spec{2, 1}(range,:,:)*exp(-2*1i*pi/6), 1))); % надо в районе 5 кгц и 2 кгц

% contourf(data.x, data.z, viz); 

p = pcolor(data.x, data.z, viz);
% p = pcolor(viz);
set(p, "LineStyle", "none", "FaceColor", "interp");
colormap jet;


%%
range=[300:350];
clf;
stpz=data.z(1,2)-data.z(1,1);
beta_vect=[-1/stpz/2:1/squeeze((data.z(1,end)-data.z(1,1))):1/stpz/2];
beta_spec = mean(mean((fft(data.spec{2, 1}(range,5:10,:), [],3)), 1), 2); 
plot((beta_vect),fftshift(squeeze(abs(beta_spec))));
title(['f= ',num2str(data.f([min(range),max(range)])'/1e3,'%.2f - %.2f'),' kHz'])
grid on
grid minor

%%
plot(data.z(5,:), mean(viz(5:8,:), 1), "LineWidth", 1); % показать (корелляционная функция)
set(gca, "FontSize", 16);
xlabel("z, mm");
ylabel("C_{xy}");
title("Correlation function");
grid on;

%% графики автоспектров (точка pos = (5,2))

z_coord = 10;

pxx1 = real(data.spec{1,1}(:, z_coord, 24));
pxx2 = real(data.spec{2,2}(:, z_coord, 24));

figure();
%subplot(2,1,1);
loglog(f, pxx1);
xlabel('Hz');
ylabel('');
title('CTA');
grid on;
xlim([0 fs/2]);

hold on;

%subplot(2,1,2);
loglog(f, pxx2);
xlabel('Hz');
ylabel('');
title('Surface sensor');
grid on;
xlim([0 fs/2]);

% взаимный спектр (точка pos = (5,2))

Sxy = data.spec{1,2}(:, z_coord, 24);

pxy = abs(Sxy);
phi = angle(Sxy);

figure();
subplot(2,1,1);
loglog(f, pxy);
xlabel('Hz');
ylabel('|Gxy|');
title('амплитуда');
grid on;
xlim([1 fs/2]);

subplot(2,1,2);
plot(f, rad2deg(unwrap(phi)));
xlabel('Hz');
ylabel('градусы');
title('фаза');
grid on;
xlim([1 fs/2]);
%ylim([-180 180]);

% разобраться до конца со взаимным спектром, построить для разных точек
% и попытаться найти фазовую скорость по наклону графика

%% гифка: sweep по фазе с шагом 10 градусов
% формула: exp(-n*1i*pi/18), шаг 10° => n = 0,1,...,35

% 5-15 (250Hz) ; 100-150 (3KHz) ; 300-350 (8KHz)
range = 300:350;
gif_filename = 'phase_sweep_8khz.gif';

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
    title(sprintf('\\phi = %d°  (n = %d)', n*10, n));
    xlabel('x, mm');
    ylabel('z, mm');
    drawnow;

    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if n == 0
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end
