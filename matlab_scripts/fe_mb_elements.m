elements = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe'};
base_dir = '/data/margo/atomic_calculations/fe_elements/RSCAN_norm_val';
recip_grid = linspace(0,52.917720860738790,2001).*1.889/(2*pi);
Z = 1:length(elements);

for i = 1:length(elements)
    element_dir = fullfile(base_dir, elements{i});
    
    try
        fx_fe_plot_val_norm(element_dir);
        close all   % only reached if no error
    catch ME
        warning('Skipping %s due to error: %s', elements{i}, ME.message);
        continue;   % go to next loop iteration
    end
    
    %fe_hf = fe_lobato_tab(recip_grid, i);
    %fe_out = [recip_grid' fe_hf'];
    %save(fullfile(element_dir, 'fe_hf.xy'), 'fe_out', "-ascii");
end
%%

HF_DFT_MSE_LDA = [];

for i = 1:length(elements)
    element_dir = fullfile(base_dir, elements(i));
    fe_hf = load(fullfile(element_dir, 'fe_hf.xy'));
    fe_dft = load(fullfile(element_dir, 'fe_ae.xy'));
    rmse_vals = sqrt(mean((fe_hf(:,2) - fe_dft(:,2)).^2));
    HF_DFT_MSE_LDA = [HF_DFT_MSE_LDA , rmse_vals];

    figure(i)
    hold on
    box on
    plot(fe_hf(:,1), fe_hf(:,2), 'b--', LineWidth=2)
    plot(fe_dft(:,1), fe_dft(:,2), 'k-', LineWidth=2)
    xlim([0, 3])
    legend('f_e HF', 'f_e DFT', 'FontSize', 14)
    legend boxoff
    saveas(gcf, fullfile(element_dir, 'fe_hf_dft.png'))

end

%%
figure(19)
hold on
box on
plot(Z, HF_DFT_MSE_LDA, 'ks-')
plot(Z, HF_DFT_MSE_PBE, 'bs-')
plot(Z, HF_DFT_MSE_RSCAN, 'rs-')
legend("LDA", "PBE", "RSCAN", 'FontSize', 14)
legend boxoff
xlabel("Z", "FontSize",14)
ylabel("RMSE", "FontSize",14)

%%
base_dir_lda = "/data/margo/atomic_calculations/fe_elements/LDA"
base_dir_pbe = "/data/margo/atomic_calculations/fe_elements/PBE"
base_dir_rscan = "/data/margo/atomic_calculations/fe_elements/RSCAN" 

STD_ae = []
STD_core = []
STD_val = []

for i = 1:length(elements)

    element_dir_lda = fullfile(base_dir_lda, elements(i));
    element_dir_pbe = fullfile(base_dir_pbe, elements(i));
    element_dir_rscan = fullfile(base_dir_rscan, elements(i));


    fe_lda_ae = load(fullfile(element_dir_lda, 'fe_ae.xy'));
    fe_lda_core = load(fullfile(element_dir_lda, 'fe_core.xy'));
    fe_lda_val = load(fullfile(element_dir_lda, 'fe_val.xy'));

    fe_pbe_ae = load(fullfile(element_dir_pbe, 'fe_ae.xy'));
    fe_pbe_core = load(fullfile(element_dir_pbe, 'fe_core.xy'));
    fe_pbe_val = load(fullfile(element_dir_pbe, 'fe_val.xy'));

    fe_rscan_ae = load(fullfile(element_dir_rscan, 'fe_ae.xy'));
    fe_rscan_core = load(fullfile(element_dir_rscan, 'fe_core.xy'));
    fe_rscan_val = load(fullfile(element_dir_rscan, 'fe_val.xy'));

    ae = [fe_lda_ae(:,2), fe_pbe_ae(:,2), fe_rscan_ae(:,2)];
    core = [fe_lda_core(:,2), fe_pbe_core(:,2), fe_rscan_core(:,2)];
    val = [fe_lda_val(:,2), fe_pbe_val(:,2), fe_rscan_val(:,2)];

    std_ae_val = std(ae, 0, 2)';
    std_core_val = std(core, 0, 2)';
    std_val_val = std(val, 0, 2)';

    STD_ae = [STD_ae; std_ae_val];
    STD_core = [STD_core; std_core_val];
    STD_val = [STD_val; std_val_val];     

end 

%%
figure(1)
h = heatmap(recip_grid(1:50), 1:size(STD_val,1), STD_val(:,1:50));

% Convert x-values to scientific notation
labels = compose('%.1e', recip_grid(1:50));

% Skip every second label
labels(2:2:end) = {''};   % replaces every 2nd label with empty string

% Apply formatted labels
h.XDisplayLabels = labels;
h.Title = 'Standard deviation between LDA, PBE, and RSCAN';
h.YLabel = 'Atomic Number';
h.XLabel = 'q [A^{-1}]';
h.FontSize = 14;
%%
figure(2)
h = heatmap(recip_grid(1:50), 1:size(STD_ae,1), STD_ae(:,1:50));

% Convert x-values to scientific notation
labels = compose('%.1e', recip_grid(1:50));

% Skip every second label
labels(2:2:end) = {''};   % replaces every 2nd label with empty string

% Apply formatted labels
h.XDisplayLabels = labels;
h.Title = 'Standard deviation between LDA, PBE, and RSCAN';
h.YLabel = 'Atomic Number';
h.XLabel = 'q [A^{-1}]';
h.FontSize = 14;
%%
figure(3)
h = heatmap(recip_grid(1:50), 1:size(STD_core,1), STD_core(:,1:50));

% Convert x-values to scientific notation
labels = compose('%.1e', recip_grid(1:50));

% Skip every second label
labels(2:2:end) = {''};   % replaces every 2nd label with empty string

% Apply formatted labels
h.XDisplayLabels = labels;
h.Title = 'Standard deviation between LDA, PBE, and RSCAN';
h.YLabel = 'Atomic Number';
h.XLabel = 'q [A^{-1}]';
h.FontSize = 14;
%%
base_dir = '/home/margo/PROJP_Timstem_FILES/hBN_exp_data/STEM_Opt/hBN_calc/LDA/SCF';


b = load(fullfile(base_dir, 'phase_b.txt'));     
n = length(b);               
semb = std(b) / sqrt(n)  

nt = load(fullfile(base_dir, 'phase_n.txt'));               
semn = std(nt) / sqrt(n) 

av_e = (semb+semn)/20.

%% Recalculating the N- electron scattering factor from IUCrJ (2018). 5, 348–353 Koji Yonekura et al.  Ionic scattering factors of the atoms of biological molecules 
recip_grid = linspace(0,52.917720860738790,2001).*1.889/(2*pi);

s = recip_grid./ (2)

a = [2.192e-1 7.256e-1 1.398e+0 1.245e+0 4.381e-1]
b = [4.784e-1 3.389e+0 1.171e+1 3.604e+1 1.125e+2]

a1 = [0.0267 0.1328 0.5301 1.1020 0.4215] 
b1 = [0.0541 0.5165 2.8207 10.6297 34.3764]

a = a(:)
b = b(:)


a1 = a1(:)
b1 = b1(:)

f_e_N_m = a.' *exp(-(b/4).*(recip_grid.^2)) - 4*0.02393366*(1./recip_grid.^2)
f_e_N = a1.' *exp(-b1.*(s.^2))
fe_out = [recip_grid' f_e_N_m'];
fe_n_lobato = fe_lobato_tab(recip_grid, 7) %
 
figure(10)
hold on
plot(recip_grid, f_e_N_m, 'b--')
plot(recip_grid, f_e_N, 'b')
plot(recip_grid, fe_n_lobato, 'r')


xlim([0.05 2.5])
ylim([-5 5])
save('/data/margo/PROJP_Timstem_FILES/elements_ionised/N/fe_N_ae_1neg.xy', 'fe_out', "-ascii")

%% Comparingh the GL parametrised ionic and neutral scattering factors for N 
recip_grid = linspace(0,52.917720860738790,2001).*1.889/(2*pi);
a_0 = 0.529177210544;

a = [3.21317199252197E-15	0.149083956484541	6.39628426567761E-14]';
b = [687.63682352118	0.0797254354199858	31382.7579719734]';
c = [0.207132165928444	0.85865640421877	1.09734240105222]';
d = [0.141105846936543	1.19486684308037	4.98281424804057]';

fe_ionic = fe_GL(recip_grid, a, b, c, d) - (1/(2*pi^2*a_0)) * (1 ./ (recip_grid).^2);

a_1 = [1.54358893295244E-11	0.197499834617669	5.70401074959514E-12]';
b_1 = [1518.70457738402	9.51800370646769	3805.32360270674]';
c_1 = [0.00172949031401361	0.000562095522186712	0.00203645234159251]';
d_1 = [0.163443236131514	0.0208116981142898	0.0652633732711663]';

a_2 = [6.2203176179791E-12	1.47795024890251E-12	0.482741611698911]';
b_2 = [14671.6575657141	3905035.57030734	0.261845180393484]';
c_2 = [0.19032156181347	0.0291582055959045	0.206464298810272]';
d_2 = [1.29517363445344	0.421782725763005	3.3185944443169]';

fe_atomic = fe_GL(recip_grid, a_1, b_1, c_1, d_1) + fe_GL(recip_grid, a_2, b_2, c_2, d_2)
fe_lobato = fe_lobato_tab(recip_grid, 7) %

figure(1)
hold on
plot(recip_grid, f_e_N_m, 'k.-')
plot(recip_grid, fe_ionic, 'r')
plot(recip_grid, fe_lobato, 'b')

xlim([0.00 3.5])
ylim([-1.2 6])



%% Comparingh the GL parametrised ionic and neutral scattering factors for B and the kappa scattering factor

fe_b_lobato = fe_lobato_tab(recip_grid, 5)
fe_b_ionic = fe_GL(recip_grid, [0.0625332337819677	0.327895162002454	1.86055566747295E-15], [0.237216041646524	1.61797662071637	59.3139112548067], [0.463044672264353	0.305822727571693	0.145438675811369], [3.10411712003169	6.78745707824061	1.28195850842971]) +  (1/(2*pi^2*a_0)) * (1 ./ recip_grid.^2)

fe_kappa = fe_GL(recip_grid, [3.95464212622568E-15	8.69016764967438E-12	0.195421633214958], [3.2761164056363	6872.44282204361	4.50868449677155], [0.00346698671572955	0.0040002147327075	0.000718189782390323], [0.106595899276725	0.286271786321969	0.0311185269286683] ) + 2*fe_GL(recip_grid, [0.0964801683925392	8.77397969459654E-14	3.50305379096725E-14] , [0.134758532908817	65961.1282698477	63102.7159780487] ,[0.104097892330072	0.115192698858846	0.0154549653417483], [2.87513066795704	7.05050107331932	0.976941920462288]) +  (1/(2*pi^2*a_0)) * (1 ./ recip_grid.^2)
fe_data_b_ion_numerical = load('/data/margo/PROJP_Timstem_FILES/elements_ionised/B/fe_val.xy')
fe_b_nm = fe_data_b_ion_numerical(:,2) +  (1/(2*pi^2*a_0)) * (1 ./ fe_data_b_ion_numerical(:,1).^2);
fe_b_ionic_G = fe_G(recip_grid, [0.5278	0.6091	0.1070	0.1753	0.0565], [25.7335	8.8169	0.0376	1.9779	0.3522]/4)+  0.09573464 * (1 ./ recip_grid.^2)

%figure(2)
hold on
box on 
plot(recip_grid, fe_b_ionic_G - 0.1, 'k')
plot(recip_grid, fe_b_lobato, 'b')
%plot(recip_grid, fe_kappa, 'k')
plot(recip_grid, fe_b_nm, 'm--')
%plot(s, fe_b_test, 'c')
xlim([0.13 3])
%ylim([-1.2 6])
legend('B^+ ionic (DFT)', 'B neutral (HF)', 'B^+ \kappa formalism', 'FontSize', 14)
legend boxoff
xlabel('q')
ylabel('f_e(q)')
saveas(gcf, 'fe_B_ionic_neutral_kappa.png')


g = fe_data_b_ion_numerical(:,1);
s = g/2

fe_b_numerical = fe_data_b_ion_numerical(:,2)%+  (1/(8*pi^2*a_0))* (1 ./ s.^2);
fe_b_test = fe_G(s, [0.1298 0.5200 1.3767 0.7666 ], [0.4180 3.9214 16.6634 51.7511 ])


figure(102)
hold on
plot(s, fe_b_numerical, 'm')
plot(s, fe_b_test, 'c')
xlim([0.04 2.5])

figure(101)
hold on
plot(recip_grid, fe_b_lobato, 'r')
plot(recip_grid, fe_b_test, 'b--')
xlim([0 2.5])



%% fe_Be2+ fe_O2+ vs K-fe for Be and O

recip_grid = linspace(0,52.917720860738790,2001).*1.889/(2*pi);


fe_be_ionic_G = fe_G(recip_grid, [3.40E-03	1.03E-02	2.33E-02	3.25E-02	1.20E-02], [2.67E-02	1.62E-01	5.31E-01	1.48E+00	3.88E+00])+  0.09573464 * (2 ./ recip_grid.^2)
fe_be_kappa = fe_GL(recip_grid, [2.3715494676727E-164	3.43085718003405E-166	3.71479353218328E-165], [8.39910907	4.82920858	11.82638622],[1.91289344039486E-169	2.85665385133798E-169	1.12950436689075E-169], [6.899171826	1.809554354	23.28124175]) +  0.09573464 * (2 ./ recip_grid.^2)

fe_o_ionic_G = fe_G(recip_grid, [4.21E-02	2.10E-01	8.52E-01	1.82E+00	1.17E+00], [6.09E-02	5.59E-01	2.96E+00	1.15E+01	3.77E+01]) -  0.09573464 * (2 ./ recip_grid.^2)
fe_o_kappa = fe_GL(recip_grid, [5.75556154341272E-12, 0.200620319, 5.97780028697153E-14 ], [567.7424031	12.41232454	638.6084795], [0.001844672	3.91392298240516E-13	0.000922748], [0.087587707	23.26413117	0.02623686])+(4/3)*fe_GL(recip_grid, [3.75525727780392E-13	0.579197354	1.13997006177469E-15], [449.2784266	0.341878793	53108.84055], [0.148941702	0.022023749	0.171280036], [0.913826747	0.29258548	2.367586208])-  0.09573464 * (2 ./ recip_grid.^2)
%fe_o_kappa_095 = fe_GL(recip_grid, [5.75556154341272E-12, 0.200620319, 5.97780028697153E-14 ], [567.7424031	12.41232454	638.6084795], [0.001844672	3.91392298240516E-13	0.000922748], [0.087587707	23.26413117	0.02623686])+(4/3)*fe_GL(recip_grid, [3.75525727780392E-13	0.579197354	1.13997006177469E-15], [449.2784266	0.341878793	53108.84055], [0.148941702	0.022023749	0.171280036], [0.913826747	0.29258548	2.367586208])-  0.09573464 * (2 ./ recip_grid.^2)



figure(1)
hold on
box on
plot(recip_grid, fe_be_ionic_G, 'k')
plot(recip_grid, fe_be_kappa, 'b')
xlim([0.1, 2 ])
ylim([0, 5])
legend('Be^{2+} f_e(q) ionic', 'Be^{2+} \kappa - f_e(q) ', 'FontSize', 14)
legend boxoff
xlabel('q [Å^{-1}]', 'FontSize', 14)
ylabel('f_e(q) ', 'FontSize', 14)
saveas(gcf, 'fe_Be.png')



figure(2)
hold on
box on
plot(recip_grid, fe_o_ionic_G, 'k')
plot(recip_grid, fe_o_kappa, 'b')
%plot(recip_grid, fe_o_kappa_095, 'b--')
xlim([0.1, 2 ])
ylim([-5, 2])
legend('O^{2-} f_e(q) ionic', 'O^{2-} \kappa = 1  f_e(q) ', 'O^{2-} \kappa = 1.05  f_e(q) ', 'FontSize', 14, 'Location', 'se')
legend boxoff
xlabel('q [Å^{-1}]', 'FontSize', 14)
ylabel('f_e(q) ', 'FontSize', 14)
saveas(gcf, 'fe_O.png')


%% fx_Be2+ fx_O2+ vs K-fe for Be and O


recip_grid = linspace(0,52.917720860738790,2001).*1.889/(2*pi);

fx_be_ionic_G = fe_G(recip_grid, [6.26030 0.884900 0.799300 0.164700 0], [0.002700 0.831300 2.27580 5.11460 0]) ;
fx_be_ionic_G = fx_be_ionic_G - 6.1092 ;
fx_be_kappa = load('/data/margo/atomic_calculations/fe_elements/RSCAN/Be_core/recip_core_radial_charge.out');

fx_o_ionic = fe_G(recip_grid, [3.75040  2.84294 1.54298 1.62091 0 ], [16.5151 6.59203 0.319201 43.3486 0]);
fx_o_ionic = fx_o_ionic + 0.242060 ;
fx_o_kappa_c = load('/data/margo/atomic_calculations/fe_elements/RSCAN/O/recip_core_radial_charge.out');
fx_o_kappa_v = fe_GL(recip_grid/0.9, [2.92999933098968	-2.31935004051211	-0.898904399833499	], [0.377102460556051	1.78985707894551	237.990442663738	], [0.757117413948573	-0.919670427040065	-0.3080950412633]	,[2.40038323603449	0.462431719424496	0.188962318975981]);

fx_o_kappa = fx_o_kappa_c + (4/3)*fx_o_kappa_v ; 

figure(3)
hold on
box on
plot(recip_grid, fx_be_ionic_G, 'k')
plot(recip_grid, fx_be_kappa, 'b--')
xlim([0.0, 2 ])
ylim([0, 3])
legend('Be^{2+} f_x(q) ionic', 'Be^{2+} \kappa - f_x(q) ', 'FontSize', 14)
legend boxoff
xlabel('q [Å^{-1}]', 'FontSize', 14)
ylabel('f_x(q) ', 'FontSize', 14)
saveas(gcf, 'fx_Be.png')



figure(4)
hold on
box on
plot(recip_grid, fx_o_ionic, 'k')
plot(recip_grid, fx_o_kappa, 'b--')
xlim([0.0, 2 ])
%ylim([-5, 2])
legend('O^{2-} f_x(q) ionic', 'O^{2-} \kappa = 0.9  f_x(q) ', 'FontSize', 14, 'Location', 'ne')
legend boxoff
xlabel('q [Å^{-1}]', 'FontSize', 14)
ylabel('f_x(q) ', 'FontSize', 14)
saveas(gcf, 'fx_O.png')

xlabel('q [Å^{-1}]', 'FontSize', 14)
ylabel('f_x(q) ', 'FontSize', 14)
saveas(gcf, 'fx_O.png')


%% comparing analytical O1- vs kappa O1- vs alpha mixing O1-


recip_grid = linspace(0,52.917720860738790,2001).*1.889/(2*pi);
fe_o_neutral = fe_lobato_tab(recip_grid, 8) ;
fe_o2_ionic_G = fe_G(recip_grid, [4.21E-02	2.10E-01	8.52E-01	1.82E+00	1.17E+00], [6.09E-02	5.59E-01	2.96E+00	1.15E+01	3.77E+01]) -  0.09573464 * (2 ./ recip_grid.^2)-  0.09573464 * (1 ./ recip_grid.^2) ; 
fe_o1_ionic_G = fe_G(recip_grid, [2.236E-1 6.923E-1 1.176E+0 9.354E-1 2.821E-1], [4.372E-1 2.918E+0 9.670E+0 2.868E+1 8.489E+1]) -  0.09573464 * (1 ./ recip_grid.^2) ; 

fe_mix = 0.5*fe_o2_ionic_G+0.5*fe_o_neutral;

fe_o_kappa1 = fe_GL(recip_grid, [5.75556154341272E-12, 0.200620319, 5.97780028697153E-14 ], [567.7424031	12.41232454	638.6084795], [0.001844672	3.91392298240516E-13	0.000922748], [0.087587707	23.26413117	0.02623686])+(7/6)*fe_GL(recip_grid, [3.75525727780392E-13	0.579197354	1.13997006177469E-15], [449.2784266	0.341878793	53108.84055], [0.148941702	0.022023749	0.171280036], [0.913826747	0.29258548	2.367586208])-  0.09573464 * (1 ./ recip_grid.^2) ;

fe_o_core = fe_GL(recip_grid, [5.75556154341272E-12, 0.200620319, 5.97780028697153E-14 ], [567.7424031	12.41232454	638.6084795], [0.001844672	3.91392298240516E-13	0.000922748], [0.087587707	23.26413117	0.02623686]);

fx_o_v = fe_GL(recip_grid/0.9, [2.92999933098968	-2.31935004051211	-0.898904399833499], [0.377102460556051	1.78985707894551	237.990442663738],[0.757117413948573	-0.919670427040065	-0.3080950412633],[2.40038323603449	0.462431719424496	0.188962318975981]) ;

fe_fx_v = mott_bethe(recip_grid, (7/6)*fx_o_v, 7) ;

fe_fx_o = fe_o_core + fe_fx_v -  0.09573464 * (1 ./ recip_grid.^2) ;

figure(5)
hold on
box on
plot(recip_grid, fe_o_neutral)
plot(recip_grid, fe_o1_ionic_G)
plot(recip_grid, fe_mix)
plot(recip_grid, fe_o_kappa1)
plot(recip_grid, fe_fx_o)
xlim([0.1, 2 ])
ylim([-5, 2])
legend('O neutral','O^{1-} f_e(q) ionic', 'O^{1-} mixed ', 'O^{1-} \kappa  ', 'O^{1-} \kappa = 0.95 ', 'FontSize', 14, 'Location', 'se')
legend boxoff
xlabel('q [Å^{-1}]', 'FontSize', 14)
ylabel('f_e(q) ', 'FontSize', 14)





%% O2- fe ionic vs. kappa scaled from fx properly









