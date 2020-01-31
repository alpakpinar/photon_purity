#!/usr/bin/env python
import os
import sys
from klepto.archives import dir_archive
from template_fit import make_templates,fit_templates, plot_purity, plot_comparison
pjoin = os.path.join


def main():
    inpath = os.path.abspath(sys.argv[1])
    if not os.path.isdir(inpath):
        raise IOError(f'Input directory does not exist: {inpath}')

    # Make templates if necessary
    outdir = f'./intermediate/{os.path.basename(inpath)}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    template_file = pjoin(outdir,'templates.root')

    if not os.path.exists(template_file):
        acc = dir_archive(
                    inpath,
                    serialized=True,
                    compression=0,
                    memsize=1e3,
                    )
        make_templates(acc, fout=template_file)

    # Do the fit
    # Makes fits if necessary
    variations = ['vfine','fine','nominal','coarse','vcoarse']
    for variation in variations:
        result_file = pjoin(outdir,f'results_{variation}.pkl')
        if not os.path.exists(result_file):
            fit_templates(template_file,variation)

        # Make result plots
        plot_purity(result_file,variation)
    plot_comparison(pjoin(outdir, 'results_{variation}.pkl'),variations)

if __name__ == "__main__":
    main()