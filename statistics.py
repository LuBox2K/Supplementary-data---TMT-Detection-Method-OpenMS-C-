import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def k_s_stat(sample1,sample2):
    ksstat = np.max(np. abs(sample1 - sample2))
    return ksstat

def k_s_significant(sample1,sample2,alpha):
    ks = k_s_stat(sample1,sample2)
    # We are assuming, that sample 1 and sample 2 are the same size
    crit_val = np.sqrt(-np.log(alpha / 2) * 1/sample1.size)
    return ks > crit_val , ks, crit_val # Debugging Purposes maybe?

# Generate, given a Histogram, a Cummulative Distribution, so that we can apply the K_S_Test
def get_cum_dist(array):
    #counts,edges = hist  since numpy stores hists just as tuples, this should work (verified) -> Not needed since we go to arrays
    cum_dist = np.cumsum(array)
    return cum_dist/np.sum(array) # We want to return the normalized variant, see https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#Two-sample_Kolmogorov%E2%80%93Smirnov_test

def ks_two_hists(arr1,arr2):
    cum_dist1 = get_cum_dist(arr1)
    cum_dist2 = get_cum_dist(arr2)

    d_val = k_s_stat(cum_dist1, cum_dist2)

    n1 = np.sum(arr1)
    n2 = np.sum(arr2)
    en =(n1 * n2) / (n1+n2)  # For this one look into https://docs.scipy.org/doc/scipy/tutorial/stats/continuous_kstwobign.html

    p_val = stats.kstwobign.sf(d_val * np.sqrt(en))

    return d_val, p_val


# THIS FUNCTION ENTERS CLANKERLAND. IT IS **ONLY** USED FOR PLOTTING, AS THERE IS NO FURTHER LOGIC GOING ON HERE
def plot_tmt_results(hist_rep, hist_bg, cdf_rep, cdf_bg, d_val, p_val):
    counts_r, edges = hist_rep
    counts_b, _ = hist_bg
    x_axis = edges[1:]  # Rechte Kanten für die CDF
    bin_centers = (edges[:-1] + edges[1:]) / 2 # Zentren für das Histogramm

    # Erstelle zwei Subplots untereinander
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # --- OBEN: Das klassische Histogramm ---
    width = (edges[1] - edges[0]) * 0.4
    ax1.bar(bin_centers - width/2, counts_r, width=width, label='Reporter Peaks', color='royalblue', alpha=0.8)
    ax1.bar(bin_centers + width/2, counts_b, width=width, label='Background Peaks', color='gray', alpha=0.5)
    ax1.set_ylabel('Absolute Counts')
    ax1.set_title(f'TMT Ad-Hoc Analysis (p-value: {p_val:.2e})')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # --- UNTEN: Die kumulative Verteilung (CDF) ---
    ax2.step(x_axis, cdf_rep, where='post', label='CDF Reporter', color='royalblue', lw=2)
    ax2.step(x_axis, cdf_bg, where='post', label='CDF Background', color='gray', linestyle='--', lw=2)

    # Visualisierung der maximalen Distanz D
    max_idx = np.argmax(np.abs(cdf_rep - cdf_bg))
    ax2.annotate('', xy=(x_axis[max_idx], cdf_bg[max_idx]),
                 xytext=(x_axis[max_idx], cdf_rep[max_idx]),
                 arrowprops=dict(arrowstyle='<->', color='red', lw=2))
    ax2.text(x_axis[max_idx], (cdf_rep[max_idx] + cdf_bg[max_idx])/2, f' D={d_val:.2f}', color='red', fontweight='bold')

    ax2.set_ylabel('Cumulative Probability')
    ax2.set_xlabel('m/z delta within Bin')
    ax2.set_ylim(0, 1.05)
    ax2.legend(loc='upper left')
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()

# Since we'll have 35 channels and 35 Background Channels, a wise idea might be to generate an array of p-vals.
# Needed to connect our two codes:
# A Tuple (BG, Rep) of Lists that contain 35 Lists each
# Each List represents the counts of peaks in that intervall. E.g. BG[3] contains the count of the 3+1 = 4th Background Location.

def get_pvals(backgrounds,reporters):
    res = np.zeros(shape=(35,35))
    for i in range (35):
        tmp = []
        for j in range (35):
            #dval,pval = ks_two_hists(reporters[i],backgrounds[j])
            _, pval = stats.ks_2samp(reporters[i], backgrounds[j])
            tmp.append(pval)
        res[i] = tmp
    return res

# Now that we have a Matrix of p-values, it might generally be a good Idea as to what the interpretation might be.
# We will extract multiple information first, such as the median of the pvals of one reporter and the arithmetic mean
# If there's enough time it might also be interesting to consider multiple testung, such as Bonferroni or Bonferroni Holm.
# But first...
def get_medians(matrix):
    res = np.zeros(35)
    for i in range (35):
        # Not as efficient but def. more readable.
        res[i] = np.median(matrix[i])
    return res

def get_means(matrix):
    res = np.zeros(35)
    for i in range (35):
        # Not as efficient but def more readable.
        res[i] = np.mean(matrix[i])
    return res

# Now, ongoing with the plotting mission. ALSO CLANKERLAND
def plot_sigs(bgs, reps, alpha):
    # 1. Daten vorbereiten
    p_matrix = get_pvals(bgs, reps)
    means = get_means(p_matrix)
    medians = get_medians(p_matrix)

    # Wir verteilen 35 Reporter auf Seiten à 10 Plots
    reporters_per_page = 10
    total_reporters = 35

    for page_start in range(0, total_reporters, reporters_per_page):
        # Erstelle ein Grid mit 5 Zeilen und 2 Spalten
        fig, axes = plt.subplots(5, 2, figsize=(12, 18))
        axes = axes.flatten() # Macht aus dem 5x2 Grid eine einfache Liste

        for i in range(reporters_per_page):
            rep_idx = page_start + i
            if rep_idx >= total_reporters:
                axes[i].axis('off') # Verstecke leere Subplots auf der letzten Seite
                continue

            p_values = p_matrix[rep_idx]
            sig_count = np.sum(p_values < alpha)
            nonsig_count = len(p_values) - sig_count

            # Plot in den jeweiligen Subplot zeichnen
            ax = axes[i]
            bars = ax.bar(['Sig', 'Non-Sig'], [sig_count, nonsig_count],
                          color=['#2ecc71', '#e74c3c'], alpha=0.8)

            # Werte auf die Balken schreiben
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, height + 0.2,
                        int(height), ha='center', va='bottom')

            ax.set_title(f"Rep {rep_idx+1} | Mean: {means[rep_idx]:.3f} | Med: {medians[rep_idx]:.3f}", fontsize=10)
            ax.set_ylim(0, 42)
            ax.tick_params(axis='both', which='major', labelsize=8)

        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        fig.suptitle(f"TMT Validation Page {page_start // 10 + 1}", fontsize=16, fontweight='bold')
        plt.show()

def main():
    bin_edges = np.linspace(126.1262, 126.1292, 21)

    all_reporters = []
    all_backgrounds = []

    for i in range(35):
        bg_data = np.random.uniform(126.1262, 126.1292, size=100)
        counts, _ = np.histogram(bg_data, bins=bin_edges)
        all_backgrounds.append(counts)

        if i < 20:
            rep_data = np.random.normal(loc=126.1277, scale=0.0004, size=100)
        else:
            rep_data = np.random.uniform(126.1262, 126.1292, size=100)
        counts, _ = np.histogram(rep_data, bins=bin_edges)
        all_reporters.append(counts)


    d_stat, p_val = ks_two_hists(all_reporters[0], all_backgrounds[0])

    print(f"Test mit reinen Arrays: D={d_stat:.4f}, p={p_val:.4e}")


    plot_sigs(all_backgrounds, all_reporters, alpha=0.05)

if __name__ == "__main__":
    main()
