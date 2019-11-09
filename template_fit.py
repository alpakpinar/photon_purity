import os
import re
from collections import namedtuple

import cloudpickle as pickle
import numpy as np
import uproot
from coffea import hist
from coffea.hist.export import export1d
from matplotlib import pyplot as plt

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

    pt_ax = hist.Bin('pt','$p_{T}$ (GeV)',list(range(200,400,50)) + list(range(400,700,100))  + [1000])
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

def fit_templates(template_file):
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
    for pt_tag in set(pt_tags):
        for year in [2016, 2017, 2018]:
            x = RooRealVar("x","#sigma_{i#eta i#eta}",0,0.015,"")

            def get_hist(name):
                h = f.Get(f'{name}_{pt_tag}').Clone()
                print(name, h.Integral())
                h.Rebin(4)

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

            # print(rhf_bad.dataHist().sum())
            # print(rhf_good.dataHist().sum())
            x.setRange('acceptance', 0,  0.01015 )
            x.setRange('total', 0, 0.15)
            # print('INTEGRAL:', rhf_good.createIntegral(RooArgSet(x),RooFit.Range('acceptance')).getVal())
            # print('INTEGRAL:', rhf_bad.createIntegral(RooArgSet(x),RooFit.Range('acceptance')).getVal())
            # print('INTEGRAL:', dh_data.createIntegral(RooArgSet(x),RooFit.Range('acceptance')).getVal())

            def get_facc(rhf):
                in_acc = rhf.createIntegral(RooArgSet(x),RooFit.Range('acceptance')).getVal()
                tot    = rhf.createIntegral(RooArgSet(x),RooFit.Range('total')).getVal()
                facc  = in_acc / tot
                return facc

            facc_good = get_facc(rhf_good)
            facc_bad = get_facc(rhf_bad)
            purity_in_acc = 1 / ( 1 + (facc_bad / facc_good) *(1-purity.getVal())/purity.getVal())
            chi2 = RooChi2Var("chi2","chi2",model,dh_data);
            npar=1
            t = []
            t.append(add_text(0.15,0.5,0.7,0.8, f'Purity full range = {purity.getVal():.3f} #pm {purity.getError():.2e}'))
            t.append(add_text(0.15,0.5,0.8,0.9, f'Purity in acceptance = {purity_in_acc:.3f}'))
            # t.append(add_text(0.5,0.7,0.5,0.7,"Chi2 / NDF = {0} / {1}".format(chi2.getValV(),dh_data.NumEntries()-npar)))
            # t.append(add_text(0.15,0.4,0.5,0.7,"Chi2 / NDF = {0:.1f}".format(chi2.getValV()/dh_data.numEntries())))
            # tt=add_text(0.15,0.4,0.4,0.6, f'Purity2 = {purity2.getVal():.2f}^{{ +{purity2.getErrorHi():.2f}}}_{{ {purity2.getErrorLo():.2f}}}')
            c1.SetLogy(1)
            c1.SaveAs(pjoin(outdir,f"fit_{year}_{pt_tag}.png"))

            results.append(result(
                purity_in_acceptance = purity_in_acc,
                purity_total = purity.getVal(),
                pt = pt_tag,
                year=year
            ))

    outdir = os.path.dirname(template_file)
    with open(pjoin(outdir, 'results.pkl'),'wb') as f:
        pickle.dump(results, f)

def plot_purity(result_file):
    '''Plot photon purity as a function of pt for different years.'''
    with open(result_file, 'rb') as f:
        results = pickle.load(f)

    for year in [2017,2018]:
        x, y, y_tot = [], [], []
        for result in sorted(filter(lambda x: x.year==year, results), key=lambda x: float(x.pt.replace('pt','').split('-')[0])):
            low, high = (float(x) for x in result.pt.replace('pt','').split('-'))
            x.append(0.5*(low+high))
            y.append(100 * (1-result.purity_in_acceptance))
            y_tot.append(100 * (1-result.purity_total))

        p = plt.plot(x,y,'o-',label=f'{year}, $\sigma_{{i\eta i\eta}} < 0.01015$')

        # plt.plot(x,y_tot,'o--',label=f'{year}, $\sigma_{{i\eta i\eta}} < 0.015$', color=p[0].get_color(), fillstyle='none')
        # print(x,y,y_tot)
    plt.gca().grid(1,linestyle='--')
    plt.gca().set_ylabel("Impurity (%)")
    plt.gca().set_xlabel("Photon $p_{T}$ (GeV)")
    plt.gca().set_ylim(0,5)
    plt.legend()

    outdir = pjoin('./plots/',os.path.basename(os.path.dirname(result_file)))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    plt.gcf().savefig(pjoin(outdir,'purity.pdf'))

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
