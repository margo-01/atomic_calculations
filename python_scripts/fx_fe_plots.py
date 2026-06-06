import matplotlib.pyplot as plt


def fe_single_plot(recip_grid, fe_data, filename):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(recip_grid, fe_data, "k-", linewidth=2)
    ax.set_xlim(0, 3)
    ax.xlabel(r"$g\;[\mathrm{\AA}^{-1}]$", fontsize=14)
    ax.ylabel(r"$f_e(g)\;[\mathrm{\AA}]$", fontsize=14)

    ax.legend(fontsize=14, frameon=False)
    ax.set_box_aspect(0.7)
    fig.tight_layout()
    plt.close(fig)


def fe_ae_core_val_plot(recip_grid,fe_core, fe_val,fe_ae ):
    filename = 'fe_ae_core_val_plot.png'

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(recip_grid, fe_core, "b-", linewidth=2, label="core")
    ax.plot(recip_grid, fe_val, "r-", linewidth=2, label="val")
    ax.plot(recip_grid, fe_ae, color=(128 / 255, 0, 128 / 255), linewidth=2, label="all-electron")
    ax.xlabel(r"$g\;[\mathrm{\AA}^{-1}]$", fontsize=14)
    ax.ylabel(r"$f_e(g)\;[\mathrm{\AA}]$", fontsize=14)
    ax.set_xlim(0, 3)
    ax.legend(fontsize=14, frameon=False)
    ax.set_box_aspect(0.7)
    fig.tight_layout()
    plt.close(fig)
    
def fx_single_plot(recip_grid, fx_data, filename):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(recip_grid, fx_data, "k-", linewidth=2)
    ax.set_xlim(0, 6)
    ax.xlabel(r"$g\;[\mathrm{\AA}^{-1}]$", fontsize=14)

    ax.legend(fontsize=14, frameon=False)
    ax.set_box_aspect(0.7)
    fig.tight_layout()
    plt.close(fig)

def fx_ae_core_val_plot(recip_grid, fx_core, fx_val, fx_ae, Z):
    filename = 'fx_ae_core_val_plot.png'
    
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(recip_grid, fx_core, "b-", linewidth=2, label="core")
    ax.plot(recip_grid, fx_val, "r-", linewidth=2, label="valence")
    ax.plot(recip_grid, fx_ae, color=(128 / 255, 0, 128 / 255), linewidth=2, label="all-electron")
    ax.xlabel(r"$g\;[\mathrm{\AA}^{-1}]$", fontsize=14)
    ax.ylabel(r"$f_x(g)$", fontsize=14)
    ax.set_xlim(0, 6)
    ax.set_ylim(0, Z + 1)
    ax.legend(fontsize=14, frameon=False)
    ax.set_box_aspect(0.7)
    fig.tight_layout()
    plt.close(fig)























