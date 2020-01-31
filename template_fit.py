import os
import re
from collections import namedtuple
import array
import cloudpickle as pickle
import numpy as np
import uproot
from coffea import hist
from coffea.hist.export import export1d
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import textwrap
import ROOT as r
r.gSystem.Load('libRooFit')
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from ROOT import (RooArgList, RooArgSet, RooChi2Var, RooDataHist, RooDataSet,
                  RooFit, RooHistFunc, RooHistPdf, RooRealSumPdf, RooRealVar)

pjoin = os.path.join


r.gROOT.SetBatch(r.kTRUE)
result = namedtuple('result',
                            [
                            'purity_in_acceptance',
                            'purity_total',
                            'year',
                            'pt'
                            ])

def fitfun(x, a, b, c):
    return a * np.exp(-b * x) + c

def make_templates(acc, fout):
    '''Reads coffea histograms and converts to ROOT templates.'''
    # Load inputs
    acc.load('sieie')
    acc.load('nevents')
    acc.load('sumw')

    # Scaling
    h = acc['sieie']
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    pt_ax = hist.Bin('pt','$p_{T}$ (GeV)',list(range(200,300,50)) + list(range(300,700,100))  + [1000])
    h = h.rebin('pt', pt_ax)
    h_iso = h.project('cat','medium_nosieie')
    h_noniso = h.project('cat','medium_nosieie_invertiso')

    # Make templates
    templates = {}
    for year in [2016, 2017, 2018]:
        mc = re.compile(f'(GJet).*HT.*{year}')
        data = re.compile(f'(EGamma).*{year}.*')
        templates[f'{year}_good'] = h_iso[mc].integrate('dataset')
        templates[f'{year}_bad']  = h_noniso[data].integrate('dataset')
        templates[f'{year}_data'] = h_iso[data].integrate('dataset')


    print(templates)
    # Save output
    f = uproot.recreate(fout)
    for name, histo in templates.items():
        edges = histo.axis('pt').edges()

        for i in range(len(edges)-1):
            low = edges[i]
            high = edges[i+1]

            th1 = export1d(histo.integrate('pt', slice(low, high)))
            f[f'{name}_pt{low:.0f}-{high:.0f}'] = th1

def pretty_title(pt_tag, year):
    '''Makes a nice plot title for this bin and year'''
    m = re.match(f'pt(\d+)-(\d+)', str(pt_tag))
    lo, hi = m.groups()

    return f'{year}: Photon p_{{T}}: {lo} - {hi} GeV'

def fit_templates(template_file, variation):
    '''Uses the given good/bad/data templates to perform purity fits.'''

    outdir = pjoin('./plots/', os.path.basename(os.path.dirname(template_file)))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    f = r.TFile(template_file)

    pt_tags = []
    for key in f.GetListOfKeys():
        print(key.GetName())
        m = re.match('.*(pt\d+-\d+)',key.GetName())
        if not m:
            continue
        pt_tags.append(m.groups()[0])

    results = []
    limits = {
        "nominal" : (0.004,0.02),
        "coarse" : (0.004,0.02),
        "vcoarse" : (0.004,0.02),
        "fine" : (0.004,0.02),
        "vfine" : (0.0058,0.02),
    }
    for pt_tag in set(pt_tags):
        for year in [2016, 2017, 2018]:

            x = RooRealVar("x","#sigma_{i#eta i#eta}",limits[variation][0],limits[variation][1],"")

            def get_hist(name):
                h = f.Get(f'{name}_{pt_tag}').Clone()
                print(name, h.Integral())

                if variation == 'nominal':
                    h.Rebin(5)
                elif variation == 'fine':
                    h.Rebin(2)
                elif variation == 'vfine':
                    h.Rebin(1)
                elif variation == 'coarse':
                    h.Rebin(10)
                elif variation == 'vcoarse':
                    h.Rebin(25)
                elif variation == 'twobin':
                    newbins = array.array("d",[h.GetBinLowEdge(1),0.01,h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())])
                    print("NEWBINS", newbins)
                    h = h.Rebin(len(newbins)-1, h.GetName()+"_rebin",newbins)
                else:
                    raise ValueError(f"Unknown variation: {variation}")


                nbins = h.GetNbinsX()
                overflow = h.GetBinContent(nbins+1)
                doverflow = h.GetBinError(nbins+1)
                lastbin = h.GetBinContent(nbins)
                dlastbin = h.GetBinError(nbins)

                h.SetBinContent(nbins, overflow + doverflow)
                h.SetBinError(nbins, np.hypot(doverflow, dlastbin))

                h.Scale(1./h.Integral())
                return h



            dh_data = RooDataHist("data", "data", RooArgList(x), get_hist(f'{year}_data'),1) ;
            dh_bad  = RooDataHist("bad", "bad", RooArgList(x), get_hist(f'{year}_bad'),1) ;
            dh_good = RooDataHist("good", "good", RooArgList(x), get_hist(f'{year}_good'),1) ;

            rhf_good = RooHistPdf(
                                "rhf_good",
                                "rhf_good",
                                RooArgSet(x),
                                dh_good
                                );
            rhf_bad  = RooHistPdf(
                                "rhf_bad",
                                "rhf_bad",
                                RooArgSet(x),
                                dh_bad
                                );

            purity = RooRealVar('purity', 'purity', 0.95,0.7,1)
            # purity2 = RooRealVar('purity2', 'purity2', 1,0.,5)

            model = RooRealSumPdf(
                                'model',
                                'model',
                                RooArgList(rhf_good, rhf_bad),
                                RooArgList(purity))
            model.fitTo(dh_data,RooFit.Save(1), RooFit.SumW2Error(r.kTRUE));

            c1 = r.TCanvas( 'c1','c1', 200, 10, 700, 500 )

            frame = x.frame(RooFit.Title(pretty_title(pt_tag, year)));
            dh_data.plotOn(frame)
            model.plotOn(frame, RooFit.LineColor(r.kGray+1))
            model.plotOn(frame, RooFit.Components('rhf_bad'),RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed));
            model.plotOn(frame, RooFit.Components('rhf_good'),RooFit.LineStyle(r.kDashed));

            c1.cd()
            frame.Draw()
            frame.SetMinimum(1e-7)

            # print(rhf_bad.dataHist().sum())
            # print(rhf_good.dataHist().sum())
            x.setRange('acceptance', 0,  0.01 )
            x.setRange('total', 0, 0.02)
            # print('INTEGRAL:', rhf_good.createIntegral(RooArgSet(x),RooFit.Range('acceptance')).getVal())
            # print('INTEGRAL:', rhf_bad.createIntegral(RooArgSet(x),RooFit.Range('acceptance')).getVal())
            # print('INTEGRAL:', dh_data.createIntegral(RooArgSet(x),RooFit.Range('acceptance')).getVal())

            def get_facc(rhf):
                '''Ratio of events in acceptance to all events'''
                in_acc = rhf.createIntegral(RooArgSet(x),RooFit.Range('acceptance')).getVal()
                tot    = rhf.createIntegral(RooArgSet(x),RooFit.Range('total')).getVal()
                facc  = in_acc / tot
                return facc

            facc_good = get_facc(rhf_good)
            facc_bad = get_facc(rhf_bad)
            R = facc_bad / facc_good
            purity_in_acc = 1 / ( 1 + R *(1-purity.getVal())/purity.getVal())
            purity_in_acc_err = R * purity.getError() / (-purity.getVal() * R + R + purity.getVal())
            t = []
            t.append(add_text(0.5,0.8,0.7,0.8, f'Purity full range = ({100*purity.getVal():.1f} #pm {100*purity.getError():.1f}) %'))
            t.append(add_text(0.5,0.8,0.8,0.9, f'Purity in acceptance = {100*purity_in_acc:.1f} %'))
            c1.SetLogy(1)
            c1.SaveAs(pjoin(outdir,f"fit_{year}_{pt_tag}_{variation}.png"))
            c1.SaveAs(pjoin(outdir,f"fit_{year}_{pt_tag}_{variation}.pdf"))
            results.append(result(
                purity_in_acceptance = (purity_in_acc,purity_in_acc_err),
                purity_total = purity.getVal(),
                pt = pt_tag,
                year=year
            ))

    outdir = os.path.dirname(template_file)
    with open(pjoin(outdir, f'results_{variation}.pkl'),'wb') as f:
        pickle.dump(results, f)

def plot_purity(result_file, variation):
    '''Plot photon purity as a function of pt for different years.'''
    with open(result_file, 'rb') as f:
        results = pickle.load(f)

    fig = plt.gcf()
    fig.clf()
    ax = plt.gca()
    for year in [2017,2018]:
        x, y, dy = [], [], []
        for result in sorted(filter(lambda x: x.year==year, results), key=lambda x: float(x.pt.replace('pt','').split('-')[0])):
            low, high = (float(x) for x in result.pt.replace('pt','').split('-'))
            x.append(0.5*(low+high))
            y.append(100 * (1-result.purity_in_acceptance[0]))
            dy.append(100 * result.purity_in_acceptance[1])

        print(x, y, dy)
        p = ax.errorbar(x,y,dy,fmt='o-',label=f'{year}, $\sigma_{{i\eta i\eta}} < 0.01$')

        # plt.plot(x,y_tot,'o--',label=f'{year}, $\sigma_{{i\eta i\eta}} < 0.015$', color=p[0].get_color(), fillstyle='none')
        # print(x,y,y_tot)
    ax.grid(1,linestyle='--')
    ax.set_ylabel("Impurity (%)")
    ax.set_xlabel("Photon $p_{T}$ (GeV)")
    ax.set_ylim(0,8)
    ax.legend()

    outdir = pjoin('./plots/',os.path.basename(os.path.dirname(result_file)))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fig.savefig(pjoin(outdir,f'purity_{variation}.pdf'))
    fig.savefig(pjoin(outdir,f'purity_{variation}.png'))


colors = {
    "nominal" : "#08306b",
    "fine" : "#2171b5",
    "vfine" : "#6baed6",
    "band" : "#d9d9d9",
    "fit" : "red"
}
colors['vcoarse'] = colors['vfine']
colors['coarse'] = colors['fine']
def plot_comparison(result_file_template, variations):
    '''Plot comparison of variations per year.'''

    for year in [2017,2018]:
        fig = plt.gcf()
        fig.clf()
        ax = plt.gca()
        x, y, dy = {}, {}, {}
        for variation in variations:
            with open(result_file_template.format(variation=variation), 'rb') as f:
                results = pickle.load(f)
            x[variation], y[variation], dy[variation] = [], [], []
            for result in sorted(filter(lambda x: x.year==year, results), key=lambda x: float(x.pt.replace('pt','').split('-')[0])):
                low, high = (float(x) for x in result.pt.replace('pt','').split('-'))
                x[variation].append(0.5*(low+high))
                y[variation].append(100 * (1-result.purity_in_acceptance[0]))
                dy[variation].append(100 * result.purity_in_acceptance[1])

            x[variation] = np.array(x[variation])
            y[variation] = np.array(y[variation])
            dy[variation] = np.array(dy[variation])

            if variation == 'nominal':
                fmt = 'o'
            elif 'fine' in variation:
                fmt = 's'
            elif 'coarse' in variation:
                fmt = 'x'
            p = ax.errorbar(x[variation],y[variation],dy[variation],fmt=fmt,label=f'{year}, {variation}', color=colors[variation],
            linestyle='-' if variation == 'nominal' else '--',
            zorder = 10 if variation=='nominal' else 9)

            # plt.plot(x,y_tot,'o--',label=f'{year}, $\sigma_{{i\eta i\eta}} < 0.015$', color=p[0].get_color(), fillstyle='none')
            # print(x,y,y_tot)
        ax.grid(1,linestyle='--')
        ax.set_ylabel("Impurity (%)")
        ax.set_xlabel("Photon $p_{T}$ (GeV)")
        ax.set_ylim(0,6)

        # Fit
        pars, _ = curve_fit(fitfun, x["nominal"], y["nominal"], sigma=dy["nominal"], p0=[3,1e-3,1])
        xinterp = np.linspace(min(x["nominal"]),max(x["nominal"]), 1000)

        # Plot
        unc = (1.25,1.25)
        ax.fill_between(xinterp, fitfun(xinterp, *pars) / unc[0], fitfun(xinterp, *pars) * unc[1], label=f'nominal + lnN {1/unc[0]:.2f} / {unc[1]:.2f}',color=colors['band'],zorder=-2)

        ax.plot(xinterp, fitfun(xinterp, *pars),color=colors["fit"], linestyle="-",zorder=11,linewidth=2, label="Fit to nominal")

        ax.legend(title="Binning variations")

        ax.text(0.05, 0.75,
                textwrap.dedent(
                    f"""
                    $\\bf{{f(x) = a\\times exp(-b \\times x) + c}}$
                    a = {pars[0]:.2f}
                    b = {pars[1]:.2e} / GeV
                    c = {pars[2]:.2f}
                    """),
                fontsize=10,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
        # Save plots
        outdir = pjoin('./plots/',os.path.basename(os.path.dirname(result_file_template)))
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for extension in ['pdf','png']:
            fig.savefig(pjoin(outdir,f'purity_variations_{year}.{extension}'))

def add_text(x1, x2, y1, y2, TEXT, color=r.kBlack, alignment=22, angle = 0, argument="NDC", size = None):
   T = r.TPaveText(x1,y1,x2,y2, argument);
   T.SetFillColor(0);
   T.SetFillStyle(0);
   T.SetLineColor(0);
   T.SetTextAlign(alignment);
   T.SetTextColor(color);

   if (not isinstance(TEXT, str)):
      for this_text in TEXT:
         text = T.AddText(this_text);
         text.SetTextAngle(angle);
         text.SetTextAlign(alignment);
   else:
      text = T.AddText(TEXT);
      text.SetTextAngle(angle);
      text.SetTextAlign(alignment);
   T.SetTextFont(42);
   if(size):
      T.SetTextSize(size);
   T.Draw("same");
   T.SetBorderSize(0);
   return T
