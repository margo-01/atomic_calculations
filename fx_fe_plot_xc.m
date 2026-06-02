function  f = fx_fe_plot_xc(dir,Z, Z_core, Z_val)
% provide the directory as '/my/dir'
%

recip_grid = linspace(0,52.917720860738790, 2001).*1.889/(2*pi);
fx_ae_lda = load(fullfile(dir, 'LDA/recip_iam_radial_charge.out'));
fx_core_lda = load(fullfile(dir, 'LDA/recip_core_radial_charge.out'));
fx_val_lda = load(fullfile(dir, 'LDA/recip_val_radial_charge.out'));

fx_ae_pbe = load(fullfile(dir, 'PBE/recip_iam_radial_charge.out'));
fx_core_pbe = load(fullfile(dir, 'PBE/recip_core_radial_charge.out'));
fx_val_pbe = load(fullfile(dir, 'PBE/recip_val_radial_charge.out'));

fx_ae_rscan = load(fullfile(dir, 'RSCAN/recip_iam_radial_charge.out'));
fx_core_rscan = load(fullfile(dir, 'RSCAN/recip_core_radial_charge.out'));
fx_val_rscan = load(fullfile(dir, 'RSCAN/recip_val_radial_charge.out'));


figure(1)

box on
hold on
plot(recip_grid, fx_core_lda, 'b-', 'LineWidth', 2)
plot(recip_grid, fx_val_lda, 'r-', 'LineWidth', 2)
plot(recip_grid, fx_ae_lda, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

plot(recip_grid, fx_core_pbe, 'b--', 'LineWidth', 2)
plot(recip_grid, fx_val_pbe, 'r--', 'LineWidth', 2)
plot(recip_grid, fx_ae_pbe,'--', 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)


plot(recip_grid, fx_core_rscan, 'b:', 'LineWidth', 2)
plot(recip_grid, fx_val_rscan, 'r:', 'LineWidth', 2)
plot(recip_grid, fx_ae_rscan,':', 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)


legend('LDA core', 'LDA valence', 'LDA all-electron', 'PBE core', 'PBE valence', 'PBE all-electron', 'RSCAN core', 'RSCAN valence', 'RSCAN all-electron','FontSize', 14)
legend boxoff
ylim([0, (Z+1)])
xlim([0, 6])
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{x}(q) ', 'FontSize', 14)

saveas(gcf, fullfile(dir, 'fx_xc_plot'))

fe_ae_lda = mott_bethe(recip_grid, fx_ae_lda, Z);
fe_core_lda = mott_bethe(recip_grid, fx_core_lda, Z_core);
fe_val_lda = mott_bethe(recip_grid, fx_val_lda, Z_val);

fe_ae_pbe = mott_bethe(recip_grid, fx_ae_pbe, Z);
fe_core_pbe = mott_bethe(recip_grid, fx_core_pbe, Z_core);
fe_val_pbe = mott_bethe(recip_grid, fx_val_pbe, Z_val);

fe_ae_rscan = mott_bethe(recip_grid, fx_ae_rscan, Z);
fe_core_rscan = mott_bethe(recip_grid, fx_core_rscan, Z_core);
fe_val_rscan = mott_bethe(recip_grid, fx_val_rscan, Z_val);

figure(2)

box on
hold on
plot(recip_grid, fe_core_lda, 'b-', 'LineWidth', 2)
plot(recip_grid, fe_val_lda, 'r-', 'LineWidth', 2)
plot(recip_grid, fe_ae_lda, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

plot(recip_grid, fe_core_pbe, 'b--', 'LineWidth', 2)
plot(recip_grid, fe_val_pbe, 'r--', 'LineWidth', 2)
plot(recip_grid, fe_ae_pbe,'--', 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)


plot(recip_grid, fe_core_rscan, 'b:', 'LineWidth', 2)
plot(recip_grid, fe_val_rscan, 'r:', 'LineWidth', 2)
plot(recip_grid, fe_ae_rscan,':', 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)


legend('LDA core', 'LDA valence', 'LDA all-electron', 'PBE core', 'PBE valence', 'PBE all-electron', 'RSCAN core', 'RSCAN valence', 'RSCAN all-electron','FontSize', 14)
legend boxoff

xlim([0, 3])
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{e}(q) ', 'FontSize', 14)


saveas(gcf, fullfile(dir, 'fe_xc_plot'))

fe_ae_out_lda = [recip_grid' fe_ae_lda'];
fe_val_out_lda = [recip_grid' fe_val_lda'];
fe_core_out_lda = [recip_grid' fe_core_lda'];

save(fullfile(dir,'fe_ae_lda.xy'), 'fe_ae_out_lda',"-ascii")
save(fullfile(dir,'fe_val_lda.xy'), 'fe_val_out_lda',"-ascii")
save(fullfile(dir,'fe_core_lda.xy'), 'fe_core_out_lda',"-ascii")


fe_ae_out_pbe = [recip_grid' fe_ae_pbe'];
fe_val_out_pbe = [recip_grid' fe_val_pbe'];
fe_core_out_pbe = [recip_grid' fe_core_pbe'];

save(fullfile(dir,'fe_ae_pbe.xy'), 'fe_ae_out_pbe',"-ascii")
save(fullfile(dir,'fe_val_pbe.xy'), 'fe_val_out_pbe',"-ascii")
save(fullfile(dir,'fe_core_pbe.xy'), 'fe_core_out_pbe',"-ascii")


fe_ae_out_rscan = [recip_grid' fe_ae_rscan'];
fe_val_out_rscan = [recip_grid' fe_val_rscan'];
fe_core_out_rscan = [recip_grid' fe_core_rscan'];

save(fullfile(dir,'fe_ae_pbe.xy'), 'fe_ae_out_rscan',"-ascii")
save(fullfile(dir,'fe_val_pbe.xy'), 'fe_val_out_rscan',"-ascii")
save(fullfile(dir,'fe_core_pbe.xy'), 'fe_core_out_rscan',"-ascii")
end
