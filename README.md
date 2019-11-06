
## Prerequisite

* python 3
* ROOT + RooFit
* bucoffea

## Instructions

The purity calculation is triggered by a single script:

```bash
./fit_photon_purity.py ./path/to/input
```

Three steps are executed:

1. The coffea input histograms are converted into template ROOT histograms. The templates are split into the data, and two fit templates for "good" and "fake" photons. The templates are created separately for different pt ranges and data taking years. The templates are saved into a ROOT file in the `intermediate` directory.

2. The templates are read from the output file of step 1. A template fit is performed to estimate the purity. Plots for each pt bin and year showing the fit result are saved in the `plots` folder. The resulting purities from each bin are saved into `results.pkl`.

3. The `results.pkl` file is read and the impurity is plotted as a function of the photon pt.

**NOTE: Steps 1 and 2 are only executed if their respective output files do not exist yet**.
Therefore, if you want to re-run all of the steps, make sure to delete their respective output files first!