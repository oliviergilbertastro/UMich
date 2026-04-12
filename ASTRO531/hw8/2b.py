import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
from tqdm import tqdm
import mesa_web as m



hist_05_Msun = m.read_history(r"ASTRO531/hw8/MESA-Web_Job_04052666820/trimmed_history.data", as_table=True)
hist_5_Msun = m.read_history(r"ASTRO531/hw8/MESA-Web_Job_04052666822/trimmed_history.data", as_table=True)
hist_sun = m.read_history(r"ASTRO531/hw2/MESA-Web_Job_01252661968/trimmed_history.data", as_table=True)

def cut_pre_main_sequence(hist):
    """Only keep the history starting from the main sequence"""
    Lh = 10**hist["log_LH"]
    L  = 10**hist["log_L"]
    # ZAMS = when H burning dominates total luminosity
    zams_index = np.argmax(Lh / L > 0.99)
    return hist[zams_index:]

def first_contiguous_block(mask):
    """Return only the first contiguous True block of a boolean mask."""
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return mask  # empty
    
    start = idx[0]
    
    # find where continuity breaks
    diffs = np.diff(idx)
    breaks = np.where(diffs > 1)[0]
    
    if len(breaks) > 0:
        end = idx[breaks[0]]
    else:
        end = idx[-1]
    
    new_mask = np.zeros_like(mask, dtype=bool)
    new_mask[start:end+1] = True
    
    return new_mask


def get_evolution_tracks(hist, end_tracks=None):
    """Returns the sections of the history file corresponding to the different evolutionary phases."""
    Xc = np.array(hist["center_h1"])
    Yc = np.array(hist["center_he4"])
    Cc = np.array(hist["center_c12"])
    Oc = np.array(hist["center_o16"])
    
    dX = np.gradient(Xc) # H change
    dY = np.gradient(Yc) # He change
    dC = np.gradient(Cc) # C change
    
    N = len(Xc)
    
    # 1. Main Sequence (core H burning)
    ms_mask = Xc > 1e-3
    tams = np.where(ms_mask)[0][-1]

    # 2. RGB definition (core H exhausted:)
    rgb_candidates = (
        (Xc < 1e-4)                  # H exhausted
    )

    # 3. HB (He ignition and C increasing)
    hb_candidates = (dY < -1e-6) & (dC > 1e-10) & (Yc > 1e-2)
    hb_candidates[:tams] = False
    
    if np.any(hb_candidates):
        hb_start = np.where(hb_candidates)[0][0]
        has_hb = True
    else:
        has_hb = False
    
    rgb_candidates[:tams] = False
    if has_hb:
        rgb_candidates[hb_start:] = False
    rgb_mask = first_contiguous_block(rgb_candidates)
    
    # 4. HB (core He burning)
    hb_mask = np.zeros(N, dtype=bool)
    if has_hb:
        # HB continues while He decreasing and C increasing
        hb_full = (dY < -1e-6) & (dC > 1e-10)
        hb_full[:hb_start] = False
        
        # stop when He exhausted
        hb_end_candidates = (Yc < 1e-3)
        hb_end_candidates[:hb_start] = False
        
        if np.any(hb_end_candidates):
            hb_end = np.where(hb_end_candidates)[0][0]
            hb_mask[hb_start:hb_end] = True
        else:
            hb_mask[hb_start:] = True
    
    # 5. AGB (post-He burning)
    # He exhausted + C/O core present
    agb_mask = np.zeros(N, dtype=bool)
    agb_candidates = (Yc < 1e-3) & ((Cc + Oc) > 1e-3)
    
    if np.any(agb_candidates):
        agb_start = np.where(agb_candidates)[0][0]
        
        # enforce ordering
        if has_hb:
            hb_end = np.where(hb_mask)[0][-1]
            agb_start = max(agb_start, hb_end + 1)
        agb_mask[agb_start:] = True
    
    # Make sure there's no overlap between the conditions
    rgb_mask = rgb_mask & (~ms_mask)
    hb_mask  = hb_mask  & (~ms_mask) & (~rgb_mask)
    agb_mask = agb_mask & (~ms_mask) & (~rgb_mask) & (~hb_mask)
    if end_tracks is not None:
        end_mask = np.zeros(N, dtype=bool)
        end_mask[end_tracks:] = True
        rgb_mask = rgb_mask & (~end_mask)
        hb_mask = hb_mask & (~end_mask)
        agb_mask = agb_mask & (~end_mask)


    def phase_times(hist, ms, rgb, hb, agb):
        age = np.array(hist["star_age"])
        
        def dt(mask):
            if np.sum(mask) < 2:
                return 0.0
            return np.log10(age[mask][-1] - age[mask][0])
        
        t_ms  = dt(ms)
        t_rgb = dt(rgb)
        t_hb  = dt(hb)
        t_agb = dt(agb)
        
        t_total = age[-1] - age[0]
        
        return {
            "MS": t_ms,
            "RGB": t_rgb,
            "HB": t_hb,
            "AGB": t_agb,
            "Total": t_total
        }
    print(phase_times(hist, ms_mask, rgb_mask, hb_mask, agb_mask))
    return hist[ms_mask], hist[rgb_mask], hist[hb_mask], hist[agb_mask]

ms_05_Msun = cut_pre_main_sequence(hist_05_Msun)
ms_5_Msun = cut_pre_main_sequence(hist_5_Msun)
ms_sun = cut_pre_main_sequence(hist_sun)

color_dict = {r"$1M_\odot$":"black",r"$0.5M_\odot$":"blue",r"$5M_\odot$":"red"}
ls_dict = {r"$1M_\odot$":"-",r"$0.5M_\odot$":"--",r"$5M_\odot$":":"}
hist_dict = {r"$1M_\odot$":ms_sun,r"$0.5M_\odot$":ms_05_Msun,r"$5M_\odot$":ms_5_Msun}
end_dict = {r"$1M_\odot$":12000,r"$0.5M_\odot$":4900,r"$5M_\odot$":800}

plt.figure(figsize=(18,6))
for i, mass in enumerate([r"$0.5M_\odot$",r"$1M_\odot$",r"$5M_\odot$"]):
    ax = plt.subplot(1,3,i+1)
    hist = hist_dict[mass]
    main_seq, rgb, hb, early_agb = get_evolution_tracks(hist, end_tracks=end_dict[mass])
    ax.plot(hist["log_Teff"], hist["log_L"], color="black", ls="-", lw=3, alpha=0.7)
    ax.plot(rgb["log_Teff"], rgb["log_L"], color="orange", ls="solid", lw=5, alpha=1, label="RGB")
    ax.plot(early_agb["log_Teff"], early_agb["log_L"], color="red", ls="solid", lw=5, alpha=1, label="Early AGB")
    ax.plot(hb["log_Teff"], hb["log_L"], color="blue", ls="solid", lw=5, alpha=1, label="HB")
    ax.plot(main_seq["log_Teff"], main_seq["log_L"], color="green", ls="solid", lw=3, alpha=1, label="MS")
    ax.set_xlabel(r"$\log T_\mathrm{eff}$ [K]")
    ax.set_ylabel(r"$\log L$ [$L_\odot$]")
    ax.legend(fontsize=15, loc=("lower left" if i!=2 else "center left"))
    ax.set_title(mass, fontsize=16)
    #ax.set_ylim(-2,5.25)
    #ax.set_xlim(3.4,5.1)
    ax.invert_xaxis()
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/2b.pdf")
plt.show()






