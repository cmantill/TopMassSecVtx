#! /usr/bin/env python                                                                                                                                                                                                            
import ROOT
from array import array
import optparse
import numpy as np
import os,sys
from UserCode.TopMassSecVtx.PlotUtils import *

"""
Create the 250 pseudoexperiments, calculate the bias and fill the pull distributions  
"""
def pseudo(opt):

    # Open ROOT file with nominal sample
    fIn_nominal = ROOT.TFile.Open('distributions/nominal/MC8TeV_TTJets_MSDecays_172v5.root')
    gen_nominal = fIn_nominal.Get('ptll_gen')
    gen_nominal.SetDirectory(0)
    rec_nominal = fIn_nominal.Get('ptll_rec')
    rec_nominal.SetDirectory(0)
    mig = fIn_nominal.Get('ptll_migration')
    mig.SetDirectory(0)
    fIn_nominal.Close()

    bins_gen = [15.13,30,41.76,51.57,60.23,70.22,79.93,89.23,104.48,136.92]

    hbias=ROOT.TH2D('bias','',len(bins_gen)-1,array('d',bins_gen),10,-0.1,0.1)

    # Get number of events from properly scaled input files for nominal mass variations to estimate the actual statistical error? 
    nevtsSeed = rec_nominal.Integral()
    nevtsToGen = ROOT.gRandom.Poisson(nevtsSeed)
    nevtsToGen = int(nevtsToGen)

    pseudorec, pseudogen = None, None
    pseudorec = rec_nominal.Clone()
    pseudogen = gen_nominal.Clone()

    # Loop to generate toy experiments
    for num in range (0,250):
        pseudorec.Reset('ICE')
        pseudorec.FillRandom(rec_nominal, nevtsToGen)

        # Unfold PE
        tunfold=ROOT.TUnfoldSys(mig, ROOT.TUnfold.kHistMapOutputHoriz, ROOT.TUnfold.kRegModeCurvature)
        tunfold.SetInput(pseudorec)
        tunfold.SetBias(gen_nominal)
        scaleBias = 1.0
        tau_value = 0.0006
        tunfold.DoUnfold(tau_value, pseudorec, scaleBias)
        rec_unfolded=gen_nominal.Clone('rec_unfolded')
        rec_unfolded.Reset('ICE')
        tunfold.GetOutput(rec_unfolded)

        # Fill bias 2D histogram: 
        for xbin in xrange(1,gen_nominal.GetNbinsX()+1):
            hbias.Fill(gen_nominal.GetXaxis().GetBinCenter(xbin), float(1-rec_unfolded.GetBinContent(xbin)/gen_nominal.GetBinContent(xbin)))
        
        fOut=ROOT.TFile.Open('distributions/testpe.root','RECREATE')
        hbias.Write()
        print 'Bias histogram saved in %s' % fOut.GetName()
        fOut.Close()
       
        print "Number of pseudoexperiment:", num
   
    # Graph to fill to show the expected bias and uncertainty
    biasGraph=ROOT.TGraphErrors()
    biasGraph.SetMarkerStyle(20)
    biasGraph.SetFillStyle(0)
    
    mins = 99999
    maxs = -1
    biasAvg = 0
    biasSigma = 0
    biasAvgUnc = 0
    biasSigmaUnc = 0

    c=TCanvas('c','c')

    # Loop over 2D hbias to project bin by bin the results
    for xbin in xrange(0,hbias.GetNbinsX()):
        firstxbin = 0
        lastxbin = -1

        if xbin > 0: 
            firstxbin=xbin
            lastxbin=xbin
        histproj= hbias.ProjectionX('_py', firstxbin, lastxbin)

        if(histproj.Integral()<10): continue
        histproj.Fit("gaus","MQ+")
        gaus = ROOT.TF1(histproj.GetFunction("gaus"))
        avg = gaus.GetParameter(1)
        sigma = gaus.GetParameter(2)

        gaus.SetParLimits(1,avg-sigma,avg+sigma)
        histproj.Fit(gaus,"MQ+")
        avg   = gaus.GetParameter(1)
        sigma = gaus.GetParameter(2)

        if xbin == 0:
            biasAvg = avg
            biasSigma = sigma
            biasAvgUnc = gaus.GetParError(1)
            biasSigmaUnc = gaus.GetParError(2)

        else:
            low = int(histproj.GetXaxis().GetBinLowEdge(xbin))
            if low < mins: mins = low
            if low > maxs: maxs = low
        
            np = biasGraph.GetN()
            biasGraph.SetPoint(np,low,avg)
            biasGraph.SetPointError(np,0,sigma)

        histproj.Delete()
        
        
    ROOT.gStyle.SetOptStat(0)

    biasGr=ROOT.TGraphErrors()
    biasGr.SetFillColor(ROOT.kGray)
    biasGr.SetFillStyle(1001)
    biasGr.SetLineColor(ROOT.kBlue)
    biasGr.SetMarkerColor(ROOT.kBlue)
    biasGr.SetPoint(0,mins,biasAvg-biasSigma)
    biasGr.SetPoint(1,mins,biasAvg+biasSigma)
    biasGr.SetPoint(2,maxs,biasAvg+biasSigma)
    biasGr.SetPoint(3,maxs,biasAvg-biasSigma)
    biasGr.SetPoint(4,maxs,biasAvg-biasSigma)

    c.Clear()
    biasGr.Draw("af")
    biasGr.GetYaxis().SetRangeUser(biasGr.GetYaxis().GetXmin()*1.5,biasGr.GetYaxis().GetXmax()*1.5)
    biasGr.GetXaxis().SetTitle("P_{t}(l^{+}l^{-})")
    biasGr.GetYaxis().SetTitle("<Bias>")
    line = ROOT.TLine(mins,biasAvg,maxs,biasAvg)
    line.Draw()
    biasGraph.Draw("p")
    c.SaveAs("distributions/testPE.root")
    c.SaveAs("distributions/testPE.png")
    """
    #c1=TCanvas('c1','c1')
    #biasGr.GetYaxis().SetTitleSize(0.047)
    #biasGr.GetYaxis().SetLabelSize(0.047)
    #biasGr.GetXaxis().SetTitle("")
    #biasGr.GetYaxis().SetTitle("")
    #biasGr.Draw("p")
    #biasGr.SaveAs("distributions/PE/bias.root")
    #c1.SaveAs("distributions/PE/bias.png")
    """
def main():
  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('-v', '--var',
                          dest='var',
                          default='ptpos',
                          help='Variable to unfold (note requires var_rec,var_gen,var_migration plots stored in the ROOT file [default: %default]')
  parser.add_option('-o', '--output',
                          dest='outDir',
                          default='distributions/PE/rootfiles/',
                          help='Output directory [default: %default]')

  (opt, args) = parser.parse_args()

  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  ROOT.gROOT.SetBatch(True)
  setTDRStyle()
  ROOT.gSystem.Load("libUserCodeTopMassSecVtx")
  ROOT.AutoLibraryLoader.enable()
  os.system('mkdir -p %s' % opt.outDir)

  print 80*'-'
  pseudo(opt)
  print 80*'-'

if __name__ == "__main__":
        sys.exit(main())
