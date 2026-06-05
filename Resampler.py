import pyopenms as oms
import sys


# this is the example C++ code for the SavitzkyGolayFilter
"""
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResampler.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  auto file_dta = OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_SavitzkyGolayFilter.dta");

  // A DTA file always has exactly one Spectrum, so we get that
  MSSpectrum spectrum;
  // Load the dta file into the spectrum
  FileHandler().loadSpectrum(file_dta, spectrum);

  LinearResampler lr;
  Param param_lr;
  param_lr.setValue("spacing", 0.01);
  lr.setParameters(param_lr);
  lr.raster(spectrum);

  SavitzkyGolayFilter sg;
  Param param_sg;
  param_sg.setValue("frame_length", 21);
  param_sg.setValue("polynomial_order", 3);
  sg.setParameters(param_sg);
  sg.filter(spectrum);

  return 0;
} //end of main
"""


### args[1:] must consist of valid filepaths
def main():
    ### configuring the resampler
    lr = oms.LinearResampler()
    param = oms.Param()
    param.setValue("spacing", 0.005)
    lr.setParameters(param)

    ## iterate input paths
    for file in sys.argv[1:]:
        print(f'loading: {file}', end='')
        exp = oms.MSExperiment()
        try:
            oms.FileHandler().loadExperiment(file, exp)
            print(col( ' SUCCESS', g=255))
        except:
            print(col(' FAILED', r=255))
            exit()

        ### resample spectra
        for spectrum in exp:
            lr.raster(spectrum)

        ### storing the resampled experiment
        print(f"storing: {file[:-5]}resampled.mzML", end='')
        oms.MzMLFile().store(f"{file[:-5]}resampled.mzML", exp)
        print(col(' SUCCESS', g=255))


### makes text colorful
def col(inp, r=0, g=0, b=0):
    return f'\033[38;2;{r};{g};{b}m{inp}\033[0m'


### main hook
if __name__ == '__main__':
    main()

