function  f = fx_fe_norm_plot(dir)
% provide the directory as '/my/dir'

%
a_0 = 0.529177210544;


recip_grid = linspace(0,52.917720860738790, 2001).*1.889/(2*pi);
fx_ae = load(fullfile(dir, 'recip_iam_radial_charge.out'));
fx_core = load(fullfile(dir, 'recip_core_radial_charge.out'));
fx_val = load(fullfile(dir, 'recip_val_radial_charge.out'));
rho_ae = load(fullfile(dir,'rho_ae.out'));
rho_core = load(fullfile(dir,'rho_core.out'));
rho_val = load(fullfile(dir,'rho_val.out'));


Z = round(trapz(rho_ae(:,1), rho_ae(:,2)));

Z_core = round(trapz(rho_core(:,1), rho_core(:,2)));
disp(Z_core)
Z_val = round(trapz(rho_val(:,1), rho_val(:,2)));
disp(Z_val)

fx_val_norm = fx_val./Z_val
rho_val_norm = rho_val/Z_val


r2_mean_ae = (trapz(rho_ae(:,1), rho_ae(:,1).^2 .* rho_ae(:,2)))./(trapz(rho_ae(:,1), rho_ae(:,2)));
r2_mean_core = (trapz(rho_core(:,1), rho_core(:,1).^2 .* rho_core(:,2)))./(trapz(rho_core(:,1), rho_core(:,2)));
r2_mean_val = (trapz(rho_val(:,1), rho_val(:,1).^2 .* rho_val(:,2)))./(trapz(rho_val(:,1), rho_val(:,2)));

%r2_mean_val_norm = (trapz(rho_val_norm(:,1), rho_val_norm(:,1).^2 .* rho_val_norm(:,2)))./(trapz(rho_val_norm(:,1), rho_val_norm(:,2)));


fe_ae = mott_bethe(recip_grid, fx_ae, Z);
fe_ae(1) = (1/3*a_0)* Z * r2_mean_ae;
fe_ae(1) = interp1(recip_grid, fe_ae, recip_grid(1), 'spline', 'extrap')%(1/3*a_0)* Z * r2_mean_ae;

fe_core = mott_bethe(recip_grid, fx_core, Z_core);
fe_core(1) = (1/3*a_0)* Z_core * r2_mean_core;
fe_core(1) = interp1(recip_grid, fe_core, recip_grid(1), 'spline', 'extrap')%(1/3*a_0)* Z_core * r2_mean_core;

fe_val = mott_bethe(recip_grid, fx_val, Z_val);
fe_val(1) = (1/3*a_0)* Z_val * r2_mean_val;
fe_val(1) = interp1(recip_grid, fe_val, recip_grid(1),'spline', 'extrap' )%(1/3*a_0)* Z_val * r2_mean_val;


fe_val_norm = mott_bethe(recip_grid, fx_val_norm, 1);
fe_val_norm(1) = fe_val_norm(2);
fe_val_norm(1) = interp1(recip_grid, fe_val_norm, recip_grid(1),'spline', 'extrap' )%(1/3*a_0)* Z_val * r2_mean_val;

figure(1)

box on
hold on
plot(recip_grid, fe_core, 'b-', 'LineWidth', 2)
plot(recip_grid, fe_val, 'r-', 'LineWidth', 2)
plot(recip_grid,Z_val*fe_val_norm,  'b-.', 'LineWidth', 2)
plot(recip_grid, fe_ae, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

xlim([0, 3])

legend('core', 'val', 'sum valence', 'all-electron', 'FontSize', 14)
legend boxoff

fe_val_norm_out = [recip_grid' fe_val_norm'];

save(fullfile(dir,'fe_val_norm.xy'), 'fe_val_norm_out',"-ascii")
end

