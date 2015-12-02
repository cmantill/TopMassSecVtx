#! /usr/bin/env python
import ROOT
import utils
from array import array
import optparse
import numpy as np
import os,sys
from UserCode.TopMassSecVtx.storeTools_cff import fillFromStore
from makeSVLMassHistos import getMassTrees
from runDileptonUnfolding import *
from UserCode.TopMassSecVtx.PlotUtils import *


"""
Shows the unfolded result and saves histograms to file
"""
def showResult(data,signal,var,outDir,fname):

    outDir= os.path.join(outDir, 'plots')

    #nominal comparison
    varPlot = Plot('%s_%s'%(var,fname),False)
    varPlot.add(signal,'signal_gen',ROOT.kGray,False)
    varPlot.add(data,  'signal_unf',  1,         True) 
    varPlot.appendTo('%s/%s_%s.root' %(outDir,var,fname))
    
    #normalized results
    dataPDF=data.Clone('%s_pdf'%data.GetName())
    dataPDF.Scale(1./dataPDF.Integral('width'))
    signalPDF=signal.Clone('%s_pdf'%signal.GetName())
    signalPDF.Scale(1./signal.Integral('width'))
    varPDFPlot = Plot('%sPDF_%s'%(var,fname),False)
    varPDFPlot.add(signalPDF,'signal_gen',ROOT.kGray,False)
    varPDFPlot.add(dataPDF,  'signal_unf',  1,         True)
    varPlot.appendTo('%s/%s_%s.root' %(outDir,var,fname))
    
    #free memory
    varPlot.reset()
    varPDFPlot.reset()

"""
Unfolding cycle for pseudoexperiments
"""
def unfoldVariable(opt):

    histos={}
    
    #open ROOT file
    fIn=ROOT.TFile.Open('unfoldResults/'+opt.var+'/nominal/plots/plotter.root')

    #reconstructed level histograms
    recList       = getPathToObjects(fIn.Get('%s_rec'%(opt.var)))
    for p in recList:
        h=fIn.Get(p)
        hname=h.GetName()
        if 'TTJets' in hname:
            if not 'signal' in histos : 
                histos['signal']=h.Clone('signal')
            else : 
                histos['signal'].Add(h)
 
    #generator level histograms 
    genList       = getPathToObjects(fIn.Get('%s_gen'%(opt.var)))
    for p in genList:
        h=fIn.Get(p)
        hname=h.GetName()
        if 'TTJets' in hname:
            if not 'signal_gen' in histos : 
                histos['signal_gen']=h.Clone('signal_gen')
            else : 
                histos['signal_gen'].Add(h)

    #detach objects from file and close file
    for h in histos: histos[h].SetDirectory(0)
    fIn.Close()

    nevtsSeed = histos['signal'].Integral()
    #print nevtsSeed
    nevtsToGen = ROOT.gRandom.Poisson(nevtsSeed)
    #print nevtsToGen
    nevtsToGen = int(nevtsToGen)
    
    pseudoMC = None
    pseudoMC = histos['signal'].Clone('signal')

    hgen = None
    hgen = histos['signal_gen'].Clone('signal_gen')

    Pull=ROOT.TH1D("Pull","Pull",200,-5.,5.)

    # Generating the pseudoexperiments
    for num in range (0,10000):
      pseudoMC.Reset('ICE')
      pseudoMC.FillRandom(histos['signal'],nevtsToGen)

      #dump migration matrix in a root file or use that histogram as migration matrix
      filenamem = 'unfoldResults/'+opt.var+'/nominal/plots/migration_'+opt.var+'_TTJets.root'

      # Get migration matrix from 172v5 samples
      fInp=ROOT.TFile.Open(filenamem)       
      hmist=ROOT.TH2D() 
      hmist=fInp.Get("migration")
      hmist.SetDirectory(0)
      fInp.Close()

      # Use migration matrix read from filenamem
      tunfold=ROOT.TUnfoldSys(hmist, ROOT.TUnfold.kHistMapOutputHoriz, ROOT.TUnfold.kRegModeCurvature)

      tunfold.SetInput(pseudoMC)

      tunfold.SetBias(histos['signal_gen'])

      scaleBias = 1.0

      tau = 1e-4

      #regularization parameter
      tunfold.DoUnfold(tau, pseudoMC, scaleBias)

      #get the unfolded distribution
      data_unfolded=histos['signal_gen'].Clone('data_unfolded')
      data_unfolded.Reset('ICE')
      tunfold.GetOutput(data_unfolded)
      
      for xbin in xrange(1,data_unfolded.GetXaxis().GetNbins()+1):
            binWidth=data_unfolded.GetXaxis().GetBinWidth(xbin)
            for h in [ data_unfolded, histos['signal_gen'] ] :
                h.SetBinContent(xbin,h.GetBinContent(xbin)/binWidth)
                h.SetBinError(xbin,h.GetBinError(xbin)/binWidth)
      fname='unfolded_wgt_%s'%num

      print fname
      showResult(data=data_unfolded,
                 signal=histos['signal_gen'],
                 var=opt.var,
                 outDir=opt.output,
                 fname=fname)
    
      histos['signal_gen'] = hgen.Clone('signal_gen')

      filename = 'unfoldResults/test/plots/ptpos_unfolded_wgt_%s.root'%(num)
      fIn=ROOT.TFile.Open(filename)
      unfolded       = getPathToObjects(fIn.Get('%s_unfolded_wgt_%s'%(opt.var,num)))
      key = 'unf'
      for p in unfolded:
         h=fIn.Get(p)
         hname=h.GetName()
         if 'data_unfolded' in hname : 
             histos['unf'] = h.Clone('data_unfolded') 

      for h in histos: histos[h].SetDirectory(0)
      fIn.Close()

      mean = histos[key].GetMean()
      meanerr = histos[key].GetMeanError()
      rms = histos[key].GetRMS()
      rmserr = histos[key].GetRMSError()

      print mean

      Mtrue = 44.419352
      
      #Fill the histograms that we are interested
      pull=(Mtrue-mean)/meanerr   
      print pull
      Pull.Fill(pull)
    
    gauss = ROOT.TF1("gauss", "gaus",-5., 5.)
    gauss.SetLineColor(ROOT.kBlue)
    gauss.SetLineWidth(3)
    gauss.SetLineStyle(1)

    Pull.Fit("gauss","EI")
    #gauss.Draw()
    mn=gauss.GetParameter(1)
    mErr=gauss.GetParError(1)

    c3=TCanvas('c3','c3')
    Pull.SetMarkerStyle(8)
    Pull.GetYaxis().SetTitleSize(0.047)
    Pull.GetYaxis().SetLabelSize(0.047)
    Pull.GetXaxis().SetTitle("Pull")
    Pull.GetYaxis().SetTitle("Pseudoexperiments")
    Pull.Draw()
    Pull.SaveAs("Pull.root")
    c3.SaveAs("Pull.pdf")


def main():
  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('-i', '--input',
                          dest='input',   
                          default='/store/cmst3/group/top/summer2015/treedir_bbbcb36/ttbar/', 
                          help='input directory with the files [default: %default]')
  parser.add_option('-r', '--root',
                          dest='root', 
                          default=None,
                          help='ROOT file with distributions to unfold [default: %default]')
  parser.add_option('-v', '--var',
                          dest='var', 
                          default='ptpos',
                          help='Variable to unfold (note requires var_rec,var_gen,var_migration plots stored in the ROOT file [default: %default]')
  parser.add_option('-o', '--output',
                          dest='output', 
                          default='unfoldResults',                                                                       
                          help='Output directory [default: %default]')

  (opt, args) = parser.parse_args()

  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  ROOT.gROOT.SetBatch(True)
  setTDRStyle()
  ROOT.gSystem.Load("libUserCodeTopMassSecVtx")
  ROOT.AutoLibraryLoader.enable()
  os.system('mkdir -p %s' % opt.output)

  print 80*'-'
  print 'Unfolding variable %s from'%(opt.var)
  unfoldVariable(opt)
  print 80*'-'

if __name__ == "__main__":
	sys.exit(main())
