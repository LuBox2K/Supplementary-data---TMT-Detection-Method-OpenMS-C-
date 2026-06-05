import sys
import numpy as np
import pyopenms as oms
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp, false_discovery_control
from statistics import *
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar



## MAGICNUMBERZZ ##

tolerance = 0.0005  ## tolerance around reporter
tolerance = 0.002921999999983882/2.0  ## min channel distance/2
#minNoDeuterium = .0065
max_mz = 137.0      ## reporter region bound
min_mz = 124.0      ## reporter region bound

intensity_threshold = 1e+07
duplex            = [0, 2]
sixplex           = [0, 1, 4, 5, 8, 9]
tenplex           = range(10)
elevenplex        = range(11)
sixteenplex       = range(16)

## TODO:
## internal calibration
## lowest 16 percent in each spectrum
## filtering for the peakpicker ## SG or even better a gaussian

thirtytwoplex     = list(range(18, 34)) + list(range(16))
thrtyfiveplex     = range(35)

channel_locations = np.array([
            ## Standard 18plex
            126.127726, ## TMTpro-126    6 32 
    
            127.124761, ## TMTpro-127N   6  32
            127.131081, ## TMTpro-127C      32
    
            128.128116, ## TMTpro-128N      32
            128.134436, ## TMTpro-128C   6  32
    
            129.131471, ## TMTpro-129N   6  32
            129.137791, ## TMTpro-129C      32
    
            130.134826, ## TMTpro-130N      32
            130.141146, ## TMTpro-130C   6  32
    
            131.138181, ## TMTpro-131N   6  32 
            131.144501, ## TMTpro-131C      32
    
            132.141536, ## TMTpro-132N      32
            132.147856, ## TMTpro-132C      32
    
            133.144891, ## TMTpro-133N      32
            133.151211, ## TMTpro-133C      32
    
            134.148246, ## TMTpro-134N      32
            134.154566, ## TMTpro-134C   ##18 / 35
    
            135.151601, ## TMTpro-135N   ##18     35
            ###########################
    
            ## New Channels used in 35plex
            127.134003, ## TMTpro-127D    All in 32 below but 135CD

            128.131038, ## TMTpro-128ND
            128.137358, ## TMTpro-128CD

            129.134393, ## TMTpro-129ND
            129.140713, ## TMTpro-129CD
    
            130.137748, ## TMTpro-130ND
            130.144068, ## TMTpro-130CD
    
            131.141103, ## TMTpro-131ND
            131.147423, ## TMTpro-131CD
    
            132.144458, ## TMTpro-132ND
            132.150778, ## TMTpro-132CD
    
            133.147813, ## TMTpro-133ND
            133.154133, ## TMTpro-133CD
    
            134.151171, ## TMTpro-134ND
            134.157491, ## TMTpro-134CD

            135.154526, ## TMTpro-135ND
            135.160846, ## TMTpro-135CD   ## 35
        ])

channel_names = [
"126",
"127N", "127C", 
"128N", "128C", 
"129N", "129C", 
"130N", "130C", 
"131N", "131C", 
"132N", "132C", 
"133N", "133C", 
"134N", "134C", 
"135N", "127D",
"128ND", "128CD",
"129ND", "129CD",
"130ND", "130CD",
"131ND", "131CD",
"132ND","132CD", 
"133ND", "133CD", 
"134ND", "134CD", 
"135ND", "135CD"]

bg_locations = np.array([125.0005    , 125.32148541, 125.64297079, 125.96445618, 126.28794147,
                         126.60942686, 126.93091225, 127.25839736, 127.57988274, 127.90136813,
                         128.23085315, 128.55233854, 128.87332395, 129.20280897, 129.52429435,
                         129.84577974, 130.17526476, 130.49675015, 130.81823553, 131.14572065+.005,
                         131.46920594, 131.79069133, 132.11217671, 132.44116176, 132.76264714,
                         133.08413253, 133.41361755, 133.73510294, 134.05658832, 134.38607334,
                         134.70755873, 135.02904412, 135.35652923, 135.67801461, 135.9995     ])

#ref_peaks = []
#for ch in channel_locations:
#    rp = oms.InternalCalibration.ReferencePeak()
#    rp.mz = ch
#    ref_peaks.append(rp)
#
#calibrator = oms.InternalCalibration()
#params = calibrator.geteDefaults()
#calibrator.setParameters(params)

def main():
    global channel_locations, bg_locations

    if not sys.argv[1:]:
        print(col('no File provided', r=255))
        exit();
    else:
        file = sys.argv[1]


    ## try loading the File
    print(f'loading {file}', end='')
    exp = oms.MSExperiment()
    try:
        oms.FileHandler().loadExperiment(file, exp)
        print(col( ' SUCCESS', g=255))
    except:
        print(col(' FAILED', r=255))
        exit()

    #trim_save_exit(exp)


    maxlevel = np.max(exp.getMSLevels())
    if maxlevel <= 1:
        print(col('this is not a TMT experiment', r=255, g=255))
        exit()

    rep_accumulator = [[] for _ in range(35)] #np.array(len(exp), dtype='object')
    bg_accumulator = [[] for _ in range(35)] #np.array(len(exp), dtype='object')
    global_mz_accumulator = []
    global_int_accumulator = []


    ## accumulate spectra
    for spec_num, spectrum in enumerate(exp):
        if spectrum.getMSLevel() != maxlevel: continue

        ## Align with Internal Calibration
        #calibrator.calibrate(spectrum, ref_peaks)
        mzs, ints = spectrum.get_peaks()

        if ints.size == 0 or mzs.size == 0:
            print(col(f'spectrum nr {spec_num}: was entirely empty', r=255), file=sys.stderr)

        # create t/f vector mask for the mz range of interest
        region_mask = (mzs >= min_mz) & (mzs <= max_mz)
        mzs = mzs[region_mask]
        ints = ints[region_mask]

        if ints.size == 0 or mzs.size == 0:
            print(col(f'spectrum nr {spec_num}: reporterregion is empty', r=255), file=sys.stderr)
            continue

        ## peaks below 16% max intensity
        low16 = (ints <= np.max(ints) * .16)
        ## aggregate peaks above 16% max intensity
        global_mz_accumulator.append(mzs[~low16])
        global_int_accumulator.append(ints[~low16])

    ## flatten
    global_mz_accumulator = np.concatenate(global_mz_accumulator)
    global_int_accumulator = np.concatenate(global_int_accumulator)

    ## add ints of same mz
    global_mz_accumulator, global_int_accumulator = deduplicate_mzs_tol(global_mz_accumulator, global_int_accumulator)

    ## clean low signal on the channels
    #for i, ion in enumerate(channel_locations):
    #    near_ion = np.abs(global_mz_accumulator - float(ion)) <= tolerance
    #    too_small = (global_int_accumulator <= intensity_threshold)
    #    global_int_accumulator = global_int_accumulator[~(near_ion & too_small)]
    #    global_mz_accumulator = global_mz_accumulator[~(near_ion & too_small)]


    quick_plot(global_mz_accumulator, global_int_accumulator, 'adjusted', oms.PeakPickerHiRes())

    ## quick_plot(global_mz_accumulator, global_int_accumulator, 'IterativePicker', oms.PeakPickerIterative())
    ## quick_plot(global_mz_accumulator, global_int_accumulator, 'ChromatogramPicker', oms.PeakPickerChromatogram()) #lol
    ## quick_plot(global_mz_accumulator, global_int_accumulator, 'IMPicker', oms.PeakPickerIM())

    generate_slices_and_print_pvals(global_mz_accumulator, global_int_accumulator)
    plt.show()

def trim_save_exit(exp): ## must be oms Experiment
    maxlevel = np.max(exp.getMSLevels())
    if maxlevel <= 1:
        print(col('this is not a TMT experiment', r=255, g=255))
        exit()

    out_exp = oms.MSExperiment()
    for spectrum in exp: ## trim the spectra
        if spectrum.getMSLevel() != maxlevel: continue

        mzs, ints = spectrum.get_peaks()

        if ints.size == 0 or mzs.size == 0:
            continue

        # create t/f vector mask for the mz range of interest
        region_mask = (mzs >= min_mz) & (mzs <= max_mz)
        mzs = mzs[region_mask]
        ints = ints[region_mask]

        if ints.size == 0 or mzs.size == 0:
            continue

        out_spec = oms.MSSpectrum()
        out_spec.setMSLevel(spectrum.getMSLevel())  # ← MS Level!
        out_spec.setRT(spectrum.getRT())  # ← Retention Time!
        out_spec.setNativeID(spectrum.getNativeID())  # ← Scan ID!
        out_spec.setPrecursors(spectrum.getPrecursors())
        out_spec.set_peaks((
            np.array(mzs, dtype='float64'),
            np.array(ints, dtype='float64')
        ))
        out_exp.addSpectrum(out_spec)


    out_exp.updateRanges()
    oms.FileHandler().storeExperiment('TMT10-Plex(PXD067886)/20250326_1_1_TRIMMED.mzML', out_exp)
    print(col('TMT10-Plex(PXD067886)/20250326_1_1_TRIMMED.mzML created', 255, 255, 0))
    exit()




def generate_slices_and_print_pvals(global_mz_accumulator, global_int_accumulator):
    # init bucket containers for background and reporter samples
    bg_slices  = [[np.array([]), np.array([])] for _ in range(35)]
    rep_slices = [[np.array([]), np.array([])] for _ in range(35)]
    # walk reporters
    for i, (ion, bg_sample) in enumerate(zip(channel_locations, bg_locations)):
        near_ion = np.abs(global_mz_accumulator - float(ion)) <= tolerance
        near_bgsample = np.abs(global_mz_accumulator - float(bg_sample)) <= tolerance

        rep_slices[i] = np.stack((global_mz_accumulator[near_ion], global_int_accumulator[near_ion]))
        bg_slices[i] = np.stack((global_mz_accumulator[near_bgsample], global_int_accumulator[near_bgsample]))

    #### save all for inspection
    #for i in range(35):
    #    #save_mzml(rep_slices[i][0], rep_slices[i][1], f'spectra/{'avg'}|{i}|ION') ## debug
    #    #save_mzml(bg_slices[i][0], bg_slices[i][1], f'spectra/{'avg'}|{i}|BG')  ## debug
    #    if len(rep_slices[i][0]) > 5:
    #        quick_plot(rep_slices[i][0], rep_slices[i][1], f'ion{i}')
    #    if len(bg_slices[i][0]) > 5:
    #        quick_plot(bg_slices[i][0], bg_slices[i][1], f'bg{i}')

    #save_mzml(np.concatenate(global_mz_accumulator), np.concatenate(global_int_accumulator), 'spectra/FULL_AVERAGE')
    #quick_plot(global_mz_accumulator, global_int_accumulator)
    #plt.show()

    p_vals = np.eye(35)
    for ch in range(35):
        for bg in range(35):
            _, p_val = ks_2samp(rep_slices[ch][1], bg_slices[bg][1])
            p_vals[ch, bg] = p_val
    nanmed = np.nanmedian(p_vals, axis=1)
    print('CALCULATED PVALUES: (↓ions | backgrouds→)\n', p_vals)
    print('SIGNIFICANT PVALUES: (↓ions | backgrouds→)\n', (p_vals < .05))
    print('SIGNIFICANT MEDIANS: (Ions→)\n', (nanmed < .05))
    print('PVAL MEDIANS with BH: (Ions→)\n', false_discovery_control(np.nan_to_num(nanmed, copy=False,nan=1.0), method='bh'))
    print('SIGNIFICANT MEDIANS with BH: (Ions→)\n', np.where(false_discovery_control(np.nan_to_num(nanmed, copy=False,nan=1.0), method='bh') < 0.05))


def col(inp, r=0, g=0, b=0):
    return f'\033[38;2;{r};{g};{b}m{inp}\033[0m'


def quick_plot(mzs, ints, title='', peakpicker=None):
    print(col(f'tolerance is set to {tolerance}', g=255))
    if len(mzs) == 0 or len(ints) == 0:
        print(col(f'{title}: empty', r=255, g=128))
        return
    if len(mzs) != len(ints):
        print(col(f'{title}: mzs/ints length mismatch {len(mzs)} vs {len(ints)}', r=255))
        return
    ##max, min = np.max(mzs), np.min(mzs) ####
    plt.figure()
    plt.stem(mzs, ints, markerfmt=' ', basefmt=' ')

    for i in range(35):
        ch = channel_locations[i]
        bg = bg_locations[i]

        plt.axvline(x = ch, color = 'red', linewidth = 0.8)
        plt.axvspan(ch-tolerance, ch+tolerance, color = 'gray', alpha = 0.3)
        plt.text(ch, 0, f'{i}|    '+channel_names[i], color="red", fontsize=10, rotation=90, va='top')

        plt.axvline(x = bg, color = 'green', linewidth = 0.8)
        plt.axvspan(bg-tolerance, bg+tolerance, color='gray', alpha=0.3)
        plt.text(bg, 0, f'bg{i}', color="green", fontsize=6, rotation=90, va='top')

    ## decalibration testing area
    plt.axvspan(channel_locations[0]-0.05, channel_locations[0]+0.05, color='yellow', alpha=0.3)

    if peakpicker:
        p_mzs, p_ints = pick_peaks(mzs, ints, peakpicker)
        if len(p_mzs) > 0:
            print(col('HIER SOLLTEST DU PEAKS SEHEN',r=255))
            # zorder=5 sorgt dafür, dass es ÜBER dem stem plot liegt
            plt.scatter(p_mzs, p_ints, color='lightgreen', s=3, zorder=5, label='Peaks')

    plt.xlabel('m/z')
    plt.ylabel('intensity')
    plt.title(title)
    plt.tight_layout()

def pick_peaks(mzs_near, ints_near, picker = oms.PeakPickerHiRes()): ## expects only reporterregion
    spec = oms.MSSpectrum()
    spec.set_peaks((mzs_near.astype(np.float64), ints_near.astype(np.float64)))
    spec.sortByPosition()

    picker = oms.PeakPickerHiRes()
    picked = oms.MSSpectrum()
    picker.pick(spec, picked)

    mzs, ints = picked.get_peaks()
    return picked.get_peaks()


def save_mzml(mzs, ints=None, filename='output'):
    if ints is None:
        ints = np.zeros(len(mzs))
    
    out_spec = oms.MSSpectrum()
    out_spec.set_peaks((
        np.array(mzs, dtype='float64'),
        np.array(ints, dtype='float64')
    ))
    out_exp = oms.MSExperiment()
    out_exp.addSpectrum(out_spec)
    out_exp.updateRanges()
    oms.FileHandler().storeExperiment(filename + '.mzML', out_exp)
    print(col(filename + '.mzML created', 255, 255, 0))


def roll_background_locations(seed=None, n=35, padding_factor=1):
    rng = np.random.default_rng(seed)
    tol = tolerance * padding_factor
    chosen = []

    for _ in range(100_000):  # max Versuche
        candidate = rng.uniform(min_mz + tol, max_mz - tol)

        ## too close to channel
        if any(np.abs(candidate - ion) <= tol * 2 for ion in channel_locations): continue
        ## too close to previously picked bg location
        if any(np.abs(candidate - bg) <= tol * 2 for bg in chosen): continue

        chosen.append(candidate)
        if len(chosen) == n:
            return np.array(chosen)

    raise ValueError(f"Nur {len(chosen)}/{n} BG-Kandidaten gefunden — Range zu klein oder tolerance zu groß?")


def min_value_distance_all_pairs(arr):
    a = np.array(arr)
    diffs = np.abs(a[:, None] - a[None, :])   ## (n x n) matrix of all differences
    np.fill_diagonal(diffs, np.inf)           ## exclude self-differences (i == j)
    return float(np.min(diffs))


def weighted_median(values, weights):
    i = np.argsort(values)
    c = np.cumsum(weights[i])
    return values[i[np.searchsorted(c, 0.5 * c[-1])]]


def find_median_decalibration(mzs, ints, tolerance = 0.05): ###   default only works for Channel 0: 126Da
    near = np.abs(mzs - float(channel_locations[0])) <= tolerance

    ## return the weighted median of the values to detect peaks.
    ## weights attenuate background
    return weighted_median(mzs[near], ints[near]) - channel_locations[0]


def find_mean_decalibration(mzs, ints, tolerance = 0.05): ###   default only works for Channel 0: 126Da
    near = np.abs(mzs - float(channel_locations[0])) <= tolerance

    ## return the weighted median of the values to detect peaks.
    ## weights attenuate background
    return np.sum(mzs[near]*ints[near]) / np.sum(ints[near]) - channel_locations[0]


def find_cubic_decalibration(mzs, ints, tolerance = 0.05):
    near = np.abs(mzs - float(channel_locations[0])) <= tolerance
    mzs_near = mzs[near]
    ints_near = ints[near]

    mzs_sorted = np.sort(mzs_near)
    ints_sorted = ints_near[np.argsort(mzs_near)]

    spline = CubicSpline(mzs_sorted, ints_sorted)

    result = minimize_scalar(
        lambda mz: -spline(mz),
        bounds = (mzs_sorted[0],mzs_sorted[-1]),
        method = 'bounded'
    )

    max_peak_mz = result.x
    shift = max_peak_mz - channel_locations[0]
    return shift


def find_peakpick_decalibration(mzs, ints, tolerance = 0.05):
    near = np.abs(mzs - float(channel_locations[0])) <= tolerance
    mzs_near, ints_near = mzs[near], ints[near]

    if (not len(mzs_near) or not len(ints_near)):
        return 0

    spec = oms.MSSpectrum()
    spec.set_peaks((mzs_near.astype(np.float64), ints_near.astype(np.float64)))
    spec.sortByPosition()

    picker = oms.PeakPickerHiRes()
    picked = oms.MSSpectrum()
    picker.pick(spec, picked)

    picked_mzs, picked_ints = picked.get_peaks()
    max_peak_mz = picked_mzs[np.argmax(picked_ints)]
    print(col(f'FOUND PEAK::: mz:{max_peak_mz}, int:{np.max(picked_ints)}', g=233))


    return max_peak_mz - channel_locations[0]


def deduplicate_mzs(mzs, ints):
    # 1. Einzigartige X-Werte finden und Mapping erstellen
    unique_mz, inverse = np.unique(mzs, return_inverse=True)
    # 2. Y-Werte basierend auf dem Mapping (inverse) summieren
    summed_ints = np.bincount(inverse, weights=ints)
    return unique_mz, summed_ints


def deduplicate_mzs_tol(mzs, ints, tol=1e-5):
    # Runde auf Toleranz-Grid statt exaktem float64 Match
    rounded = np.round(mzs / tol) * tol
    unique_mz, inverse = np.unique(rounded, return_inverse=True)
    summed_ints = np.bincount(inverse, weights=ints)
    return unique_mz, summed_ints

################################################################
################################################################

if __name__ == '__main__':
    np.set_printoptions(precision=6, suppress=True, threshold=np.inf, linewidth=400)
    #roll = roll_background_locations()
    #if np.array_equal(bg_locations, roll):
    #    print(col('rerolling only found the same solution', r=200))
    #else:
    #    print('new locations found !!\n', col(f'{roll}', g=200))
    main()
