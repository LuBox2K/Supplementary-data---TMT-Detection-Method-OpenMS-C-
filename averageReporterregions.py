import sys
import pyopenms as oms
import numpy as np
from math import floor

def main():
    ## MAGICNUMBERZZ ##
    max_mz = 150.0
    min_mz = 110.0
    bin_size = 0.005
    filelist = sys.argv[1:]

    cut_spectra = [] ## only save reporter ion region of interest

    for file in filelist:
        print(f'processing {file}', end='')
        exp = oms.MSExperiment()
        try:
            oms.FileHandler().loadExperiment(file, exp)
            print()
        except:
            print('\033[38;2;255;0;0m FAILED\033[0m')


        for spec in exp:
            if spec.getMSLevel() != 2: continue
 
            # get peaks
            mzs, ints = spec.get_peaks()

            # create t/f vector mask for the mz range
            mask = (mzs >= min_mz) & (mzs <= max_mz)
            cut_spectra.append((mzs[mask], ints[mask]))


    ## preparing bins
    bin_count = int(max_mz / bin_size) + 1
    bins = [[] for _ in range(bin_count)]

    for mzs, ints in cut_spectra:
        for m, i in zip(mzs, ints):
            bin = int(m / bin_size)
            if bin < bin_count: ## unnecessary safeguard but it might break a long computation run so....
                bins[bin].append(i)
            else:
                print(col('SOMETHING WENT VERY WRONG', 255, 0, 0))

    final_mz = []
    final_int = []

    ## calc median for each bin
    for i in range(bin_count):
        if bins[i]:
            final_mz.append(i * bin_size)
            # Median berechnen (numpy ist hier schneller)
            final_int.append(np.median(bins[i]))

    # write mzML
    out_spec = oms.MSSpectrum()
    out_spec.set_peaks((np.array(final_mz, dtype='float64'),
                        np.array(final_int, dtype='float64')))

    out_exp = oms.MSExperiment()
    out_exp.addSpectrum(out_spec)
    oms.FileHandler().storeExperiment("BACKGROUND.mzML", out_exp)
    print(col('BACKGROUND.mzML created', 255, 255, 0))



def col(str, r=0, g=0, b=0):
    return f'\033[38;2;{r};{g};{b}m{str}\033[0m'

################################################################
################################################################
if __name__ == '__main__':
    #experiment = oms.MSExperiment()
    #oms.FileHandler().loadExperiment(sys.argv[1], experiment)
    #spectrum = oms.MSSpectrum()
    #spectrum.set_peaks(experiment[0].get_peaks())
    #out = oms.MSExperiment()
    #out.addSpectrum(spectrum)
    #oms.FileHandler().storeExperiment("BACKGROUND.mzML", out)
    #input('checkmal ob das spectrum generiert wurde')

    main()

