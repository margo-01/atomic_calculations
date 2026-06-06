%Scattering factors and potentials plots

a_0 = 0.529177210544
Z = 1.00
recip_grid = linspace(0,52.917720860738790, 2001).*1.889/(2*pi)


fx_analyt_h = @(x) 1./((1+(pi^2)*(a_0^2)*(x.^2)).^2)

fx_test = fx_analyt_h(recip_grid)
figure(1)
box on
hold on
fplot(fx_analyt_h,[0 5])
%%
fe_analyt_h = @(x) (a_0/2)*((1./(1+(pi^2)*(a_0^2)*(x.^2)))+(1./(1+(pi^2)*(a_0^2)*(x.^2)).^2) ) 

fe_test = fe_analyt_h(recip_grid)
figure(1)
box on
hold on
fplot(fe_analyt_h,[0 5])
%%
recip_grid = linspace(0,52.917720860738790, 2001).*1.889/(2*pi)

fx_h_mgga = load('/home/margo/PROJP_Timstem_FILES/atomic_calculations/H/RSCAN/3/recip_iam_radial_charge.out')

fx_h_gga = load('/home/margo/PROJP_Timstem_FILES/atomic_calculations/H/PBE/2/Castep_run/recip_iam_radial_charge.out')

fx_h_lda = load('/home/margo/PROJP_Timstem_FILES/atomic_calculations/H/LDA/2/Castep_run/recip_iam_radial_charge.out')


figure(1)
box on
hold on
plot(recip_grid(1:250), fx_h_mgga(1:250))
plot(recip_grid(1:250), fx_h_gga(1:250))
plot(recip_grid(1:250), fx_h_lda(1:250))
plot(recip_grid(1:250), fx_test(1:250))

%%

fe_conv_test = mott_bethe(recip_grid, fx_test, 1)
fe_h_mgga = mott_bethe(recip_grid, fx_h_mgga, 1)

%%
a1 = 6.473848488352918e-03
a2 = -4.901925767802290e-01
a3 = 5.732841603908765e-01
a4 = -3.794033014839905e-01
a5 = 5.544264747740791e-01
b1 = 2.785198853791489e+00
b2 = 2.776204283306448e+00
b3 = 2.775385910506251e+00
b4 = 2.767593028672588e+00
b5 = 2.765118976429275e+00

a = [6.473848488352918e-03	-4.901925767802290e-01	5.732841603908765e-01	-3.794033014839905e-01	5.544264747740791e-01];
b = [2.785198853791489e+00	2.776204283306448e+00	2.775385910506251e+00	2.767593028672588e+00	2.765118976429275e+00];

test_h = fe_lobato(recip_grid, a,b)
figure(1)
box on
hold on
fe_h_param = @(x) a1*((2+b1.*x.^2)./((1+b1.*x.^2)).^2)+ a2*((2+b2.*x.^2)./((1+b2.*x.^2)).^2)+ a3*((2+b3.*x.^2)./((1+b3.*x.^2)).^2)+ a4*((2+b4.*x.^2)./((1+b4.*x.^2)).^2)+ a5*((2+b5.*x.^2)./((1+b5.*x.^2)).^2)
fe_param_test = fe_h_param(recip_grid)

plot(recip_grid, fe_param_test, 'k-')
plot(recip_grid, test_h, 'r--')
xlim([0 2])

%plot(recip_grid(1:250), fe_test(1:250), 'r--')

%%
fx_h_param = @(x) 1-2*(pi^2)*a_0.*x.^2.*fe_param_test
fx_param_test = fx_h_param(recip_grid)


%%

figure(2)
box on
hold on
plot(recip_grid(1:250), fe_conv_test(1:250))
%plot(recip_grid(1:250), fe_h_mgga(1:250))
plot(recip_grid(1:250), fe_test(1:250), '-.')

legend('test', 'mGGA', 'fe analyt')

%%
figure(1)
hold on
box on

plot(recip_grid(1:250), fx_test(1:250), 'k-', 'LineWidth', 2)
plot(recip_grid(1:250), fx_param_test(1:250), 'r:', 'LineWidth', 2)



plot(recip_grid(1:250), fe_test(1:250), 'b-', 'LineWidth', 2)
plot(recip_grid(1:250), fe_param_test(1:250), 'm--', 'LineWidth', 2)
plot(recip_grid(1:250), fe_conv_test(1:250), 'c :', 'LineWidth', 2)

legend('f_x analytical', 'f_x parametrised', 'f_e analytical', 'f_e parametrised', 'f_x analytical MB', 'FontSize', 14)
legend boxoff
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{H}(q) ', 'FontSize', 14)
saveas(gcf,'/home/margo/Desktop/param_test.png')

%%
figure(2)
box on
hold on
plot(recip_grid(1:250), fx_test(1:250), 'k-', 'LineWidth', 2)
plot(recip_grid(1:250), fx_h_mgga(1:250), 'r-', 'LineWidth', 2)
plot(recip_grid(1:250), fx_h_gga(1:250), 'm-', 'LineWidth', 2)
plot(recip_grid(1:250), fx_h_lda(1:250), 'c-', 'LineWidth', 2)

legend('f_x analytical', 'f_x mGGA', 'f_x GGA', 'f_x LDA', 'FontSize', 14)
legend boxoff
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{x,H}(q) ', 'FontSize', 14)
saveas(gcf,'/home/margo/Desktop/fx_dft_test.png')

%%

fe_h_mgga = mott_bethe(recip_grid, fx_h_mgga, 1)
fe_h_gga = mott_bethe(recip_grid, fx_h_gga, 1)
fe_h_lda = mott_bethe(recip_grid, fx_h_lda, 1)


figure(3)

box on
hold on
plot(recip_grid(1:250), fe_test(1:250), 'k-', 'LineWidth', 2)
plot(recip_grid(1:250), fe_h_mgga(1:250), 'r-', 'LineWidth', 2)
plot(recip_grid(1:250), fe_h_gga(1:250), 'm-', 'LineWidth', 2)
plot(recip_grid(1:250), fe_h_lda(1:250), 'c-', 'LineWidth', 2)

legend('f_e analytical', 'f_e mGGA', 'f_e GGA', 'f_e LDA', 'FontSize', 14)
legend boxoff
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{e,H}(q) ', 'FontSize', 14)
saveas(gcf,'/home/margo/Desktop/fe_dft_test.png')


%%
% Plotting the valence and core xrd form factors for Si LDA

fx_si_ae_lda = load('/data/margo/atomic_calculations/Si_core_val/LDA/recip_iam_radial_charge.out');

fx_si_core_lda = load('/data/margo/atomic_calculations/Si_core_val/LDA/recip_core_radial_charge.out');
fx_si_val_lda = load('/data/margo/atomic_calculations/Si_core_val/LDA/recip_val_radial_charge.out');
recip_grid = linspace(0,52.917720860738790, 2001).*1.889/(2*pi);

figure(1)
box on
hold on
plot(recip_grid, fx_si_core_lda, 'b-', 'LineWidth', 2)
plot(recip_grid, fx_si_val_lda, 'r-', 'LineWidth', 2)
plot(recip_grid, fx_si_ae_lda, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

legend('core', 'valence', 'all-electron', 'FontSize', 14)
legend boxoff
ylim([0, 15])
xlim([0, 6])
%plot(recip_grid(1:250), fx_test(1:250))
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{x,Si}(q) ', 'FontSize', 14)

% Converting valence and core scattering factors for Si from xrd to electron

fe_si_ae_lda = mott_bethe(recip_grid, fx_si_ae_lda, 6);
fe_si_core_lda = mott_bethe(recip_grid, fx_si_core_lda, 2);
fe_si_val_lda = mott_bethe(recip_grid, fx_si_val_lda, 4);

figure(2)
box on
hold on
plot(recip_grid, fe_si_core_lda, 'b-', 'LineWidth', 2)
plot(recip_grid, fe_si_val_lda, 'r-', 'LineWidth', 2)
plot(recip_grid, fe_si_ae_lda, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

xlim([0, 3])
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{e,Si}(q) ', 'FontSize', 14)

legend('core', 'val','all-electron', 'FontSize', 14)
legend boxoff

%%
%GGA

fx_si_ae_pbe = load('/data/margo/atomic_calculations/Si_core_val/PBE/recip_iam_radial_charge.out');

fx_si_core_pbe = load('/data/margo/atomic_calculations/Si_core_val/PBE/recip_core_radial_charge.out');
fx_si_val_pbe = load('/data/margo/atomic_calculations/Si_core_val/PBE/recip_val_radial_charge.out');

figure(3)
box on
hold on
plot(recip_grid, fx_si_core_pbe, 'b-', 'LineWidth', 2)
plot(recip_grid, fx_si_val_pbe, 'r-', 'LineWidth', 2)
plot(recip_grid, fx_si_ae_pbe, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

legend('core', 'valence', 'all-electron', 'FontSize', 14)
legend boxoff
ylim([0, 15])
xlim([0, 6])
%plot(recip_grid(1:250), fx_test(1:250))


% Converting valence and core scattering factors for Si from xrd to electron

fe_si_ae_pbe = mott_bethe(recip_grid, fx_si_ae_pbe, 6);
fe_si_core_pbe = mott_bethe(recip_grid, fx_si_core_pbe, 10);
fe_si_val_pbe = mott_bethe(recip_grid, fx_si_val_pbe, 4);

figure(4)
box on
hold on
plot(recip_grid, fe_si_core_pbe, 'b-', 'LineWidth', 2)
plot(recip_grid, fe_si_val_pbe, 'r-', 'LineWidth', 2)
plot(recip_grid, fe_si_ae_pbe, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

xlim([0, 3])

legend('f_{e, si} core', 'f_{e, si} val','f_{e, si} all-electron', 'FontSize', 14)
legend boxoff


%%
%mGGA

fx_si_ae_rscan = load('/data/margo/atomic_calculations/Si_core_val/RSCAN/recip_iam_radial_charge.out');

fx_si_core_rscan = load('/data/margo/atomic_calculations/Si_core_val/RSCAN/recip_core_radial_charge.out');
fx_si_val_rscan = load('/data/margo/atomic_calculations/Si_core_val/RSCAN/recip_val_radial_charge.out');

figure(5)
box on
hold on
plot(recip_grid, fx_si_core_rscan, 'b-', 'LineWidth', 2)
plot(recip_grid, fx_si_val_rscan, 'r-', 'LineWidth', 2)
plot(recip_grid, fx_si_ae_rscan, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

legend('core', 'valence', 'all-electron', 'FontSize', 14)
legend boxoff
ylim([0, 15])
xlim([0, 6])


% Converting valence and core scattering factors for Si from xrd to electron

fe_si_ae_rscan = mott_bethe(recip_grid, fx_si_ae_rscan, 6);
fe_si_core_rscan = mott_bethe(recip_grid, fx_si_core_rscan, 10);
fe_si_val_rscan = mott_bethe(recip_grid, fx_si_val_rscan, 4);

figure(6)
box on
hold on
plot(recip_grid, fe_si_core_rscan, 'b-', 'LineWidth', 2)
plot(recip_grid, fe_si_val_rscan, 'r-', 'LineWidth', 2)
plot(recip_grid, fe_si_ae_rscan, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)

xlim([0, 3])

legend('core', 'val','all-electron', 'FontSize', 14)
legend boxoff

%%

figure(7)

box on
hold on

plot(recip_grid, fx_si_ae_rscan, 'k.' ,'LineWidth', 1)
plot(recip_grid, fx_si_ae_pbe, 'r.-', 'LineWidth', 1)
plot(recip_grid, fx_si_ae_lda, 'b-', 'LineWidth', 1)

plot(recip_grid, fx_si_core_rscan, 'k.-' ,'LineWidth', 1)
plot(recip_grid, fx_si_core_pbe, 'r-', 'LineWidth', 1)
plot(recip_grid, fx_si_core_lda, 'b.', 'LineWidth', 1)

plot(recip_grid, fx_si_val_rscan, 'k-' ,'LineWidth', 1)
plot(recip_grid, fx_si_val_pbe, 'r.', 'LineWidth', 1)
plot(recip_grid, fx_si_val_lda, 'b.-', 'LineWidth', 1)



xlim([0, 3])
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{e,Si}(q) ', 'FontSize', 14)

legend('f_{e, ae} mGGA', 'f_{e, ae} GGA','f_{e, ae} LDA', 'f_{e, core} mGGA', 'f_{e, core} GGA','f_{e, core} LDA', 'f_{e, val} mGGA', 'f_{e, val} GGA','f_{e, val} LDA','FontSize', 14)
legend boxoff

%%
% Core and valence scattering factors for Si parametrised from the Zheng Wu and Zhu paper

si_ae_a = [0.746979E+0 	0.202055E+1 	0.155230E+1 	0.622582E+0 	0.413313E+0 	0.214002E+0 	0.787161E-1 	0.259177E-1]
si_ae_b = [0.806094E+2 	0.348421E+2 	0.148190E+2 	0.542240E+1 	0.183857E+1 	0.651052E+0 	0.188030E+0 	0.294949E-1]

si_core_a = [-.162443E-1 	-.200408E-2 	0.303243E-1 	0.159888E+0 	0.179614E+0 	0.932236E-1 	0.399158E-1 	0.150745E-1]
si_core_b = [0.730343E+4 	0.973802E+3 	0.593935E+1 	0.263204E+1 	0.116496E+1 	0.463095E+0 	0.145953E+0 	0.243429E-1]

si_vals_a = [0.280432E-2 	0.136895E+0 	0.434739E+0 	0.291225E+0 	0.898289E-1 	0.319712E-1 	0.979862E-2 	0.249993E-2]
si_vals_b = [0.447182E+4 	0.547416E+2 	0.245341E+2 	0.103615E+2 	0.352760E+1 	0.105564E+1 	0.268043E+0 	0.382707E-1]

si_valp_a = [0.182084E+0 	0.599125E+0 	0.495627E+0 	0.198229E+0 	0.712374E-1 	0.256984E-1 	0.834747E-2 	0.227804E-2]
si_valp_b = [0.986675E+2 	0.460214E+2 	0.205047E+2 	0.811228E+1 	0.280474E+1 	0.885983E+0 	0.236025E+0 	0.351834E-1]

fe_si_core_p = gaussian_sum(recip_grid,si_core_a, si_core_b);


fe_si_ae_p = gaussian_sum(recip_grid,si_ae_a, si_ae_b);

fe_si_valp = gaussian_sum(recip_grid,si_valp_a, si_valp_b);
fe_si_vals = gaussian_sum(recip_grid,si_vals_a, si_vals_b);

fe_si_lobato = fe_lobato(recip_grid, [2.871891426116124e+00	-2.061735011951735e+00	2.171140242044787e+00	-6.630736330588019e-02	3.010707096705137e-03], [5.084871036429896e+00	4.291781853051262e-01	3.664854341921622e-01	1.197106112969034e-01	1.439945361283975e-02]);


figure(3)
box on
hold on
plot(recip_grid, fe_si_core_p, 'b-', 'LineWidth', 2)
plot(recip_grid, (2*fe_si_valp+2*fe_si_vals), 'r-', 'LineWidth', 2)
plot(recip_grid, fe_si_ae_p, 'Color',[128, 0, 128] / 255 ,'LineWidth', 2)
legend('core', 'val','all-electron', 'FontSize', 14)
legend boxoff

xlim([0, 3])


figure(4)
box on
hold on

plot(recip_grid, fe_si_ae_p, 'b-')

plot(recip_grid, fe_si_core_p, 'b--')
plot(recip_grid, (2*fe_si_valp+2*fe_si_vals), 'b:')

plot(recip_grid, fe_si_lobato, 'k-')
plot(recip_grid, fe_si_ae_lda, 'r-' )
legend('Parametrised all-e (Zheng et al.)','Parametrised core (Zheng et al.)' ,'Parametrised valence (Zheng et al.)','Parametrised (Lobato)', 'CASTEP LDA' , 'FontSize', 14)
legend boxoff
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{e,Si}(q) ', 'FontSize', 14)
%plot(recip_grid, (fe_si_val+fe_si_core), 'k--','LineWidth', 2)
%plot(recip_grid, (2*fe_si_valp+2*fe_si_vals), 'm.')


xlim([0, 3])

%% testing the electron scattering factors that are obtained by transforming the hartree potential directly

fev_si_ae = load('/data/margo/atomic_calculations/Si_core_val/LDA/recip_ae_potential.out');

%fev_si_core = load('/data/margo/atomic_calculations/Si_core_val/LDA/recip_core_potential.out');
%fev_si_val = load('/data/margo/atomic_calculations/Si_core_val/LDA/recip_val_potential.out');

figure(4)
box on
hold on
plot(recip_grid, fev_si_ae, 'b-', 'LineWidth', 2)
%plot(recip_grid, fev_si_core, 'r-', 'LineWidth', 2)
%plot(recip_grid,  fev_si_val, 'k--','LineWidth', 2)
legend('all-electron','core','val', 'FontSize', 14)
legend boxoff
%xlim([0, 2])

%did not work out :(

%% 



rho_ae_rscan = load('/data/margo/atomic_calculations/Si_core_val/RSCAN/rho_ae.out');

rho_ae_pbe = load('/data/margo/atomic_calculations/Si_core_val/PBE/rho_ae.out');

rho_ae_lda = load('/data/margo/atomic_calculations/Si_core_val/LDA/rho_ae.out');



figure(1)
hold on
box on
plot(rho_ae_lda(:,1),rho_ae_lda(:,2), 'k.')
plot(rho_ae_pbe(:,1),rho_ae_pbe(:,2), 'b.')
plot(rho_ae_rscan(:,1),rho_ae_rscan(:,2), 'r.')
legend('LDA', 'PBE', 'RSCAN')

xlim([0, 10])

%%
%scattering factors with KH, Dirac, and SRG
recip_grid = linspace(0,52.917720860738790, 2001).*1.889/(2*pi)


fe_i_srg = load('/data/margo/atomic_calculations/Core_val_atomic/SRG/I/fe_ae_pbe.xy');
fe_i_dir = load('/data/margo/atomic_calculations/Core_val_atomic/Dirac/I/fe_ae_pbe.xy');
fe_i_kh = load('/data/margo/atomic_calculations/Core_val_atomic/KH/I_core_val/fe_ae_pbe.xy');

fe_he_srg = load('/data/margo/atomic_calculations/Core_val_atomic/SRG/He/fe_ae_pbe.xy');
fe_he_dir = load('/data/margo/atomic_calculations/Core_val_atomic/Dirac/He/fe_ae_pbe.xy');
fe_he_kh = load('/data/margo/atomic_calculations/Core_val_atomic/KH/He_core_val/fe_ae_pbe.xy');


fe_he_hf = fe_lobato(recip_grid, [3.057451160998355e+00	-6.200447791273253e+01	6.400555370846145e+01	-5.001325785427806e+00	1.517988287005264e-01], [1.089672487260788e+00	9.398387981431211e-01	9.252890343862655e-01	8.229474987086506e-01	5.773931106754022e-01] );

fe_i_hf = fe_lobato(recip_grid, [9.644006662721232e+00	-1.229244353501125e+02	1.186825648165673e+02	8.950259470255390e-02	2.612761214904750e-04], [2.726455453665441e+00	1.237234258508660e+00	1.200620368722498e+00	5.066786822851999e-02	9.554378833834973e-04] );
 
figure(1)
hold on
box on
plot(fe_i_srg(:,1), fe_i_srg(:,2), 'b-', 'LineWidth', 2)
plot(fe_i_kh(:,1), fe_i_kh(:,2), 'k-', 'LineWidth', 2)
plot(fe_i_dir(:,1), fe_i_dir(:,2), 'r-', 'LineWidth', 2)
plot(recip_grid, fe_i_hf, 'k--')
xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{e}(q) ', 'FontSize', 14)

xlim([0, 3])


legend('Schroedinger', 'Koelling-Harmon', 'Dirac', 'Relativistic HF, Lobato' , 'FontSize', 14)
legend boxoff

saveas(gcf, 'fe_rel_i.png')



figure(2)
hold on
box on
plot(fe_he_srg(:,1), fe_he_srg(:,2), 'b-', 'LineWidth', 2)
plot(fe_he_kh(:,1), fe_he_kh(:,2), 'k-', 'LineWidth', 2)
plot(fe_he_dir(:,1), fe_he_dir(:,2), 'r-', 'LineWidth', 2)
plot(recip_grid, fe_he_hf, 'k--')

xlabel('q [A^{-1}]', 'FontSize', 14)
ylabel('f_{e}(q) ', 'FontSize', 14)

xlim([0, 3])


legend('Schroedinger', 'Koelling-Harmon', 'Dirac', 'Relativistic HF, Lobato', 'FontSize', 14)
legend boxoff



saveas(gcf, 'fe_rel_he.png')


%% Calculating the fe(q) for q=0

m0 = 9.10938356e-31;    % electron mass (kg)
e  = 1.60217662e-19;    % elementary charge (C)
h  = 6.62607015e-34;    % Planck constant (J·s)
a_0 = 0.529177210544;
rho_si = load('/data/margo/atomic_calculations/Core_val_atomic/KH/Si_core_val/PBE/rho_val.out');

Z = trapz(rho_si(:,1), rho_si(:,2));

% Compute <r^2> (Eq. 6)
numerator   = trapz(rho_si(:,1), rho_si(:,1).^2 .* rho_si(:,2));
denominator = trapz(rho_si(:,1), rho_si(:,2));
r2_mean = numerator / denominator;

% Compute f_e(0) (Eq. 5)
fe0 = (1/3*a_0)* Z * r2_mean;


%% Converting the ionic scattering factors
% from this paper P. C. Schmidt and A. Weiss • X-ray Scattering Factors of Ions in Crystals

%in A^-1

q_ax = [0.0 0.02 0.04 0.05 0.06 0.08 0.10 0.12 0.14 0.15 0.16 0.18 0.2 0.25 0.3 0.35 0.40 0.45 0.50 0.60 0.70 0.80 0.90 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0];

% B3+, free ion

fx_B3 = [2.000, 1.999, 1.997, 1.995, 1.993, 1.987, 1.979, 1.970, 1.960, 1.954, 1.948, 1.934, 1.919, 1.876, 1.825, 1.767, 1.704, 1.637, 1.567, 1.420, 1.273, 1.131, 0.997, 0.875, 0.764, 0.666, 0.580, 0.505, 0.439, 0.383, 0.334, 0.292, 0.256, 0.224]

fx_N3 = [10.000 9.917 9.673 9.497 9.288 8.789 8.206 7.576 6.928 6.606 6.290 5.685 5.125 3.972 3.162 2.625 2.274 2.043 1.884 1.685 1.556 1.452 1.356 1.263 1.171 1.080 0.993 0.910 0.831 0.758 0.690 0.627 0.569 0.517];
 %[10.000, 9.133, 7.279, 6.373, 5.603, 4.532, 3.950, 3.653, 3.497, 3.445, 3.402, 3.326, 3.249, 3.016, 2.746, 2.481, 2.247, 2.054, 1.899, 1.684, 1.547, 1.448, 1.363, 1.279, 1.195, 1.110, 1.026, 0.943, 0.864, 0.789, 0.719, 0.654, 0.594, 0.539]


fe_N3 = mott_bethe_ionic(q_ax(8:end), fx_N3(8:end), 10, -3)

fe_B3 = mott_bethe_ionic(q_ax(8:end), fx_B3(8:end), 5, 3)


fe_B3_sc = mott_bethe(q_ax, fx_B3, 2)
fe_N3_sc = mott_bethe(q_ax, fx_N3, 10)


figure(2)
plot(q_ax(8:end), fe_B3_sc(8:end))

figure(3)
plot(q_ax(8:end), fe_N3_sc(8:end))

%% interpolation to fix the issues with small q behaviour

% original data (as provided)
% original data
q_ax = [0.0 0.02 0.04 0.05 0.06 0.08 0.10 0.12 0.14 0.15 0.16 0.18 0.2 0.25 0.3 0.35 0.40 0.45 0.50 0.60 0.70 0.80 0.90 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0];

fx_B3 = [2.000, 1.999, 1.997, 1.995, 1.993, 1.987, 1.979, 1.970, 1.960, 1.954, 1.948, 1.934, 1.919, 1.876, 1.825, 1.767, 1.704, 1.637, 1.567, 1.420, 1.273, 1.131, 0.997, 0.875, 0.764, 0.666, 0.580, 0.505, 0.439, 0.383, 0.334, 0.292, 0.256, 0.224];

fx_N3 = [10.000, 9.133, 7.279, 6.373, 5.603, 4.532, 3.950, 3.653, 3.497, 3.445, 3.402, 3.326, 3.249, 3.016, 2.746, 2.481, 2.247, 2.054, 1.899, 1.684, 1.547, 1.448, 1.363, 1.279, 1.195, 1.110, 1.026, 0.943, 0.864, 0.789, 0.719, 0.654, 0.594, 0.539];

% --- collapse duplicates (e.g. q=0.04) by averaging the corresponding y-values
[G, q_unique] = findgroups(q_ax);             % q_unique corresponds to group IDs
q_unique = q_unique(:);                       % make column vector

fx_B3_unique = splitapply(@mean, fx_B3, G);   % column vector (usually)
fx_N3_unique = splitapply(@mean, fx_N3, G);

% ensure column orientation
fx_B3_unique = fx_B3_unique(:);
fx_N3_unique = fx_N3_unique(:);

% --- create fine q grid (choose spacing you want)
q_fine = (0:0.001:2.0)';   % column vector, change step as desired

% --- interpolation (compute both pchip and spline variants)
fx_B3_fine_pchip  = interp1(q_unique, fx_B3_unique, q_fine, 'pchip');
fx_B3_fine_spline = interp1(q_unique, fx_B3_unique, q_fine, 'spline');

fx_N3_fine_pchip  = interp1(q_unique, fx_N3_unique, q_fine, 'pchip');
fx_N3_fine_spline = interp1(q_unique, fx_N3_unique, q_fine, 'spline');

% optional: makima if available and you want a different shape-preserving cubic
if exist('makima','file')
    fx_B3_fine_makima = makima(q_unique, fx_B3_unique, q_fine);
    fx_N3_fine_makima = makima(q_unique, fx_N3_unique, q_fine);
end

% --- quick plot to compare original vs interpolated
figure('Name','Interpolation Comparison','NumberTitle','off','Color','w');
subplot(2,1,1)
plot(q_unique, fx_B3_unique, 'o', 'DisplayName','original B3'); hold on
plot(q_fine, fx_B3_fine_pchip, '-', 'LineWidth',1, 'DisplayName','pchip interp (fine)');
plot(q_fine, fx_B3_fine_spline, '--', 'LineWidth',1, 'DisplayName','spline interp (fine)');
if exist('fx_B3_fine_makima','var')
    plot(q_fine, fx_B3_fine_makima, ':', 'LineWidth',1, 'DisplayName','makima interp (fine)');
end
xlabel('q'); ylabel('f_x B3'); legend('Location','northeast'); grid on; title('B3+ interpolation');

subplot(2,1,2)
plot(q_unique, fx_N3_unique, 'o', 'DisplayName','original N3'); hold on
plot(q_fine, fx_N3_fine_pchip, '-', 'LineWidth',1, 'DisplayName','pchip interp (fine)');
plot(q_fine, fx_N3_fine_spline, '--', 'LineWidth',1, 'DisplayName','spline interp (fine)');
if exist('fx_N3_fine_makima','var')
    plot(q_fine, fx_N3_fine_makima, ':', 'LineWidth',1, 'DisplayName','makima interp (fine)');
end
xlabel('q'); ylabel('f_x N3'); legend('Location','northeast'); grid on; title('N3 interpolation');

% --- save interpolated arrays (optional)
% save('interpolated_fx.mat','q_fine','fx_B3_fine_pchip','fx_B3_fine_spline','fx_N3_fine_pchip','fx_N3_fine_spline');

% --- Example: show first 10 values
disp('First 10 q_fine and B3 pchip values:');
disp(table(q_fine(1:10), fx_B3_fine_pchip(1:10), 'VariableNames',{'q','fx_B3_pchip'}));


fe_B3_sc_int = mott_bethe(q_fine, fx_B3_fine_pchip, 2)
fe_N3_sc_int = mott_bethe(q_fine, fx_N3_fine_pchip, 10)

