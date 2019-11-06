import re
from coffea.hist.export import export1d
from coffea import hist
import uproot
from bucoffea.plot.util import (
                                merge_extensions, 
                                scale_xs_lumi, 
                                merge_datasets
                                )


def make_templates(acc, fout):
    # Load inputs
    acc.load('sieie')
    acc.load('nevents')
    acc.load('sumw')

    # Scaling
    h = acc['sieie']
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    pt_ax = hist.Bin('pt','$p_{T}$ (GeV)',list(range(200,1000,100)) + [1200])
    h = h.rebin('pt', pt_ax)
    h_iso = h.project('cat','medium_nosieie')
    h_noniso = h.project('cat','medium_nosieie_invertiso')

    # Make templates
    templates = {}
    for year in [2017, 2018]:
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

import ROOT as r
r.gSystem.Load('libRooFit')
from ROOT import (
                  RooFit,
                  RooRealVar,
                  RooDataSet,
                  RooArgList,
                  RooDataHist,
                  RooHistFunc,
                  RooRealSumPdf,
                  RooArgSet,
                  RooChi2Var
                  )

r.gROOT.SetBatch(r.kTRUE)
def fit_templates(template_file):

    f = r.TFile(template_file)
    for key in f.GetListOfKeys():
        print(key.GetName())
    
    
    for year in [2017, 2018]:
        x = RooRealVar("x","#sigma_{i#eta i#eta}",0,0.015,"")

        def get_hist(name):
            h = f.Get(name).Clone()
            print(name, h.Integral())
            h.Rebin(2)
            # h.Scale(1./h.Integral())
            return h



        dh_data = RooDataHist("data", "data", RooArgList(x), get_hist(f'{year}_data'),1) ;
        dh_bad  = RooDataHist("bad", "bad", RooArgList(x), get_hist(f'{year}_bad'),1) ;
        dh_good = RooDataHist("good", "good", RooArgList(x), get_hist(f'{year}_good'),1) ;

        rhf_good = RooHistFunc(
                               "rhf_good",  
                               "rhf_good",
                               RooArgSet(x),
                               dh_good
                               );
        rhf_bad  = RooHistFunc(
                               "rhf_bad",
                               "rhf_bad",
                               RooArgSet(x),
                               dh_bad
                               );

        purity = RooRealVar('purity', 'purity', 1,0.,5)
        # purity2 = RooRealVar('purity2', 'purity2', 1,0.,5)

        model = RooRealSumPdf(
                            'model', 
                            'model', 
                            RooArgList(rhf_good, rhf_bad), 
                            RooArgList(purity))
        result = model.fitTo(dh_data,RooFit.Save(1), RooFit.SumW2Error(r.kTRUE));

        c1 = r.TCanvas( 'c1','c1', 200, 10, 700, 500 )
     
        frame = x.frame(RooFit.Title("title"));
        dh_data.plotOn(frame)
        model.plotOn(frame, RooFit.LineColor(r.kGray+1))
        model.plotOn(frame, RooFit.Components('rhf_bad'),RooFit.LineStyle(r.kDashed), RooFit.LineColor(r.kRed));
        model.plotOn(frame, RooFit.Components('rhf_good'),RooFit.LineStyle(r.kDashed));

        c1.cd()
        frame.Draw()


        chi2 = RooChi2Var("chi2","chi2",model,dh_data);
        npar=1
        t = []
        t.append(add_text(0.15,0.5,0.7,0.8, f'Purity = {purity.getVal():.3f} #pm {purity.getError():.2e}'))
        # t.append(add_text(0.5,0.7,0.5,0.7,"Chi2 / NDF = {0} / {1}".format(chi2.getValV(),dh_data.NumEntries()-npar)))
        t.append(add_text(0.15,0.4,0.5,0.7,"Chi2 / NDF = {0:.1f}".format(chi2.getValV()/dh_data.numEntries())))
        # tt=add_text(0.15,0.4,0.4,0.6, f'Purity2 = {purity2.getVal():.2f}^{{ +{purity2.getErrorHi():.2f}}}_{{ {purity2.getErrorLo():.2f}}}')
        c1.SetLogy(1)
        c1.SaveAs(f"fit_{year}.pdf")




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