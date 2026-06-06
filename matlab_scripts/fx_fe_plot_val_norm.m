function  f = fx_fe_plot_val_norm(dir)
% provide the directory as '/my/dir'

%
a_0 = 0.529177210544;

recip_grid = linspace(0,52.917720860738790, 2001).*1.889/(2*pi);

fx_val = load(fullfile(dir, 'recip_val_radial_charge.out'));
rho_val = load(fullfile(dir,'rho_val.out'));
Z_val = round(trapz(rho_val(:,1), rho_val(:,2)));


fx_val_norm = fx_val/Z_val ;


figure(1)
box on
hold on

plot(recip_grid, fx_val_norm, 'r-', 'LineWidth', 2)

legend( 'valence normalised', 'FontSize', 14)
legend boxoff
ylim([0, 1.2])
xlim([0, 6])
saveas(gcf, fullfile(dir, 'fx_plot_norm.png'))





r2_mean_val = (trapz(rho_val(:,1), rho_val(:,1).^2 .* rho_val(:,2)))./(trapz(rho_val(:,1), rho_val(:,2)));

r2_mean_val = r2_mean_val/Z_val ; 

fe_val_norm = mott_bethe(recip_grid, fx_val_norm, 1);
fe_val_norm(1) = (1/3*a_0)* Z_val * r2_mean_val;
fe_val_norm(1) = interp1(recip_grid, fe_val_norm, recip_grid(1),'spline', 'extrap' );%(1/3*a_0)* Z_val * r2_mean_val;

figure(2)

box on
hold on
plot(recip_grid, fe_val_norm, 'r-', 'LineWidth', 2)
xlim([0, 3])

legend('valence normalised', 'FontSize', 14)
legend boxoff
saveas(gcf, fullfile(dir, 'fe_plot_norm.png'))


fe_val_out = [recip_grid' fe_val_norm'];
fx_val_out = [recip_grid' fx_val_norm'];


save(fullfile(dir,'fe_val_norm.xy'), 'fe_val_out',"-ascii")
save(fullfile(dir,'fx_val_norm.xy'), 'fx_val_out',"-ascii")


end

