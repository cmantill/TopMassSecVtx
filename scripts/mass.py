#!/usr/bin/env python
from array import array
import ROOT
from ROOT import *
import optparse
import os,sys,numpy
import pprint

from runDileptonUnfolding import getPathToObjects
from UserCode.TopMassSecVtx.storeTools_cff import fillFromStore
from UserCode.TopMassSecVtx.PlotUtils import RatioPlot
#from runDileptonUnfolding import runDileptonUnfolding
#from runPlotter import runPlotter

"""
Get purity,stability and efficiency
"""
def getStability(root,v,m):

    histos={}
    
    #open ROOT file
    fIn=ROOT.TFile.Open(root)

    print 'analizing %s'%root

    var = str(v)

    genList       = getPathToObjects(fIn.Get('%s_gen'%(var)))
    for p in genList:
        h=fIn.Get(p)
        hname=h.GetName()
        # Signal + bkg
        if not 'gen' in histos : 
            histos['gen']=h.Clone('gen')
        else : 
            histos['gen'].Add(h)

    recList       = getPathToObjects(fIn.Get('%s_rec'%(var)))
    for p in recList:
        h=fIn.Get(p)
        hname=h.GetName()
        # Signal + bkg
        if not 'rec' in histos : 
            histos['rec']=h.Clone('rec')
        else : 
            histos['rec'].Add(h)

    migrationList = getPathToObjects(fIn.Get('%s_migration'%(var)))
    for p in migrationList:
        h=fIn.Get(p)
        hname=h.GetName()
        # Signal + bkg
        if not 'mig' in histos : 
            histos['mig']=h.Clone('mig')
        else :
            histos['mig'].Add(h)

    #pT positive lepton
    if var == 'ptpos': 
      bins_gen = [22.72,25.45,28.18,30.92,33.67,36.42,39.17,42.01,44.89,47.77,50.94,55.06,59.19,63.63,68.15,74.76,83.59,95,115.65]
      bins_rec = [20,21.36,22.73,24.09,25.45,26.82,28.18,29.55,30.92,32.29,33.67,35.04,36.41,37.79,39.17,40.57,42.01,43.45,44.89,46.33,47.77,49.21,50.93,53,55.06,57.12,59.18,61.37,63.63,65.89,68.15,70.73,74.76,78.74,83.59,88.75,95,102.85,115.65,138.46,250]
      
    title_p = ';p_{T}(l^{+}) [GeV];Purity'
    title_s = ';p_{T}(l^{+}) [GeV];Stability'
    title_e = ';p_{T}(l^{+}) [GeV];Efficiency'

    nbins_gen = histos['gen'].GetXaxis().GetNbins()
    nbins_rec = histos['mig'].GetYaxis().GetNbins()

    histos['Purity']=ROOT.TH1F('Purity',title_p,len(bins_gen)-1,array('d',bins_gen))
    histos['Stability']=ROOT.TH1F('Stability',title_s,len(bins_gen)-1,array('d',bins_gen))
    histos['Efficiency']=ROOT.TH1F('Efficiency',title_e,len(bins_gen)-1,array('d',bins_gen))

    #Ngen_tot = histos['gen'].Integral(1,nbins_gen)
    #Nrec_tot = histos['rec'].Integral(1,nbins_rec)
    #print 'gen %4.6f'%(Ngen_tot)

    # Just to check, three different ways to get the total events in hist['mig']
    #tot = histos['mig'].ProjectionY("test",1,nbins_gen,"e").Integral(1,nbins_rec)
    #tot1 = histos['mig'].Integral(1,nbins_gen,1,nbins_rec)
    #tot2 = histos['mig'].Integral()

    for i in xrange(1,len(bins_gen)):

        # Purity: # of events generated and correctly reconstructed in a given bin i relative to the # of events that are reconstructed in bin 
        # i but generated anywhere
      
        Nrec_gen = histos['mig'].Integral(i,i+1,2*i-1,2*i+1)
        Nrec = histos['mig'].Integral(1,nbins_gen,2*i-1,2*i+1)
        if Nrec != 0: 
          pur = Nrec_gen/Nrec
          histos['Purity'].SetBinContent(i,pur)

        ##  New definitions of efficiency, stability and purity
        ## Purity
        #Nrec_i = histos['rec'].GetBinContent(2*i+1) + histos['rec'].GetBinContent(2*i)
        Nrec_gen_i = histos['mig'].ProjectionY("test",i,i+1,"e").Integral(2*i,2*i+1)
        #pur2 = Nrec_gen_i/Nrec_i
        #histos['Purity'].SetBinContent(i,pur2)

        Ngen = histos['mig'].Integral(i,i+1,1,nbins_rec)
        Ngen_tot_mig = histos['mig'].Integral(1,nbins_gen,1,nbins_rec)

        #if Ngen != 0: 
        #  stab = Nrec_gen/Ngen
        #  histos['Stability'].SetBinContent(i,stab)

        #eff = Ngen/Ngen_tot_mig
        #histos['Efficiency'].SetBinContent(i,eff)

        Ngen_i = histos['gen'].GetBinContent(i)
        Ngen_rec_i = histos['mig'].ProjectionX("test2",2*i,2*i+1,"e").Integral()

        stab2 = Nrec_gen_i/Ngen_i
        histos['Stability'].SetBinContent(i,stab2)

        #eff = Nrec_gen_i/histos['rec'].Integral()
        histos['Efficiency'].SetBinContent(i, Nrec_gen_i)
         
    scale = 1.0/ histos['Efficiency'].Integral()
    histos['Efficiency'].Scale(scale)

    #norm = 1.0
    #scale = 10/histos['Efficiency'].Integral()
    #print histos['Efficiency'].GetBinContent(20)
    #print histos['Efficiency'].GetXaxis().GetBinWidth(20)
    #print histos['Efficiency'].Integral()
    #scale = histos['Efficiency'].GetXaxis().GetBinWidth(1)/ histos['Efficiency'].Integral()
    #print scale
    #histos['Efficiency'].Scale(scale)

    #dump histograms in a file\
    filename = 'quantiles/mc/%s/%s/plots.root'%(var,m)
    fOut=ROOT.TFile.Open(filename,'RECREATE')
    for h in histos: histos[h].Write()
    print 'Histograms saved in %s' % fOut.GetName()
    fOut.Close()

    # Convert histograms to TGraphErrors
    graphs ={}

    graphs['Purity'] = ROOT.TGraphErrors(histos['Purity'])
    graphs['Efficiency'] = ROOT.TGraphErrors(histos['Efficiency'])
    graphs['Stability'] = ROOT.TGraphErrors(histos['Stability'])
   
    # Save histograms in canvas

    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)


    outputName = 'quantiles/mc/'+var+'/'+m + '/'+m +'_pur_stab_eff'

    #Create a canvas for plotting your graph
    c = ROOT.TCanvas('c','c')
    gStyle.SetCanvasDefH(600);
    gStyle.SetCanvasDefW(800);
    c.SetLeftMargin(0.15);
    c.SetRightMargin(0.25)
    c.SetBottomMargin(0.25);
    p1=ROOT.TPad('p1','p1',0.,0.,1.0,1.0)
    p1.Draw()

    #Create CMS labels
    tlat1 = TLatex()
    tlat1.SetNDC()
    tlat1.SetTextFont(42)
    tlat1.SetTextSize(0.04)
    tlat1.SetTextAlign(31)
    tlat1.DrawLatex(0.16, 0.91,'#bf{CMS}')

    tlat2 = TLatex()
    tlat2.SetNDC()
    tlat2.SetTextFont(42)
    tlat2.SetTextSize(0.04)
    tlat2.SetTextAlign(31)
    tlat2.DrawLatex(0.30, 0.91, '#it{Simulation}')
    #tlat2.DrawLatex(0.38, 0.91, '#it{Simulation}')
    #tlat2.DrawLatex(0.81, 0.91, '19.7 fb^{-1}')
    #tlat2.DrawLatex(0.90, 0.91,'(8 TeV)')
    tlat2.DrawLatex(0.90, 0.91,'8 TeV')

    tleg = ROOT.TLegend(0.70, 0.75, .89, 0.89)
    tleg.SetBorderSize(0)
    tleg.SetFillColor(0)
    tleg.SetFillStyle(0)
    tleg.SetShadowColor(0)
    tleg.SetTextFont(43)
    tleg.SetTextSize(20)

    #pT positive lepton
    if var == 'ptpos': title = 'p_{T}(l^{+}) [GeV]'

    #pT charged-lepton pair
    if var == 'ptll': title = 'p_{T}(l^{+}l^{-}) [GeV]'

    #M charged-lepton pair
    if var == 'mll': title = 'M(l^{+}l^{-}) [GeV]'
      
    #Scalar sum of E
    if var == 'EposEm': title = 'E(l^{+})+E(l^{-}) [GeV]'
     
    #Scalar sum of Pt
    if var == 'ptposptm': title = 'p_{T}(l^{+})+p_{T}(l^{-}) [GeV]'

    p1.cd()
    graphs['Purity'].SetTitle("")
    graphs['Purity'].SetLineColor(kBlue-3)
    graphs['Purity'].SetLineStyle(1)
    graphs['Purity'].SetLineWidth(1)
    graphs['Purity'].SetMarkerStyle(8)
    graphs['Purity'].SetMarkerColor(kBlue-3)
    graphs['Purity'].SetMarkerSize(0.9)
    graphs['Purity'].SetMaximum(1.1)
    graphs['Purity'].SetMinimum(0.0)
    graphs['Purity'].GetXaxis().SetLabelSize(0.037)
    graphs['Purity'].GetYaxis().SetTitleOffset(1.2)
    graphs['Purity'].GetYaxis().SetTitleSize(0.037)
    graphs['Purity'].GetYaxis().SetLabelSize(0.037)
    graphs['Purity'].GetXaxis().SetTitle(title)
    graphs['Purity'].GetYaxis().SetTitle()
    for i in xrange(0,len(bins_gen)-1):
      errx = graphs['Purity'].GetErrorX(i)
      graphs['Purity'].SetPointError(i,errx,0)
    graphs['Purity'].Draw("Ap")
    tleg.AddEntry(graphs['Purity'], 'Purity', 'pl')

    graphs['Efficiency'].SetMarkerStyle(23)
    graphs['Efficiency'].SetMarkerColor(kGreen-3)
    graphs['Efficiency'].SetMarkerSize(0.9)
    graphs['Efficiency'].SetLineColor(kGreen-3)    
    graphs['Efficiency'].SetLineWidth(1)
    for i in xrange(0,len(bins_gen)-1):
      errx = graphs['Efficiency'].GetErrorX(i)
      graphs['Efficiency'].SetPointError(i,errx,0)
    graphs['Efficiency'].Draw("Samep")
    tleg.AddEntry(graphs['Efficiency'], 'Efficiency', 'pl')

    graphs['Stability'].SetMarkerStyle(21)
    graphs['Stability'].SetMarkerColor(kRed)
    graphs['Stability'].SetMarkerSize(0.9)
    graphs['Stability'].SetLineColor(kRed)    
    graphs['Stability'].SetLineWidth(1)
    for i in xrange(0,len(bins_gen)-1):
      errx = graphs['Stability'].GetErrorX(i)
      graphs['Stability'].SetPointError(i,errx,0)
    graphs['Stability'].Draw("Samep")
    tleg.AddEntry(graphs['Stability'], 'Stability', 'pl')

    tleg.Draw()

    for ext in ['png','pdf']:
      c.SaveAs('%s.%s'%(outputName,ext))
    del c


"""
Get first three Mellin moments from each of the distributions
"""
def getMoments(root,v,moments,m):

    histos={}
    
    #open ROOT file
    fIn=ROOT.TFile.Open(root)

    var = str(v)
    #unfolded histograms generated level
    if moments == 'unf':
      unfolded       = getPathToObjects(fIn.Get('%s_unfolded'%var))
      key = 'unf'
      for p in unfolded:
         h=fIn.Get(p)
         hname=h.GetName()
         if 'data_unfolded' in hname : 
             histos['unf'] = h.Clone('data_unfolded')   
 
    if moments == 'gen':
      genList       = getPathToObjects(fIn.Get('%s_gen'%(var)))
      key = 'gen'
      for p in genList:
          h=fIn.Get(p)
          hname=h.GetName()
          if 'TTJets' in hname:
              if not 'signal_gen' in histos : 
                  histos['gen']=h.Clone('gen')
              else : 
                  histos['gen'].Add(h)

    if moments == 'rec':
      recList       = getPathToObjects(fIn.Get('%s_rec'%(var)))
      key = 'signal'
      for p in recList:
          h=fIn.Get(p)
          hname=h.GetName()
          if 'Data' in hname : 
              histos['data']=h.Clone('data')
          elif 'TTJets' in hname:
              if not 'signal' in histos : 
                  histos['signal']=h.Clone('signal')
              else : 
                  histos['signal'].Add(h)
          else:
              if not 'bkg' in histos : histos['bkg']=h.Clone('bkg')
              else : histos['bkg'].Add(h)

    #detach objects from file and close file
    for h in histos: histos[h].SetDirectory(0)
    fIn.Close()

    mean = histos[key].GetMean()
    meanerr = histos[key].GetMeanError()
    rms = histos[key].GetRMS()
    rmserr = histos[key].GetRMSError()

    u0 = 1
    u1 = mean
    u1err = meanerr
    u2 = mean**2 + rms**2
    u2err = sqrt((2*rms*rmserr)**2+(2*mean*meanerr)**2)

    #Create a text file with the three different moments for the corresponding mass
    outfileName = 'quantiles/mc/'+var+'/moments_' + var + '_' + moments +'.txt'
    m = float(m)
    mass = m + 0.5
    outfile = open(outfileName,'a')
    outfile.write('%4.6f %4.6f %4.6f %4.6f %4.6f %4.6f\n' %(mass,u0,u1,u1err,u2,u2err))
    outfile.close()

from UserCode.TopMassSecVtx.PlotUtils import *

"""
Compare distributions from different masses
"""

def compareMass(compare,v):

    var = str(v)
      
    #Open root file containing histograms
    masspoints = [166,169,171,173,175,178]

    outDir = 'quantiles/mc/%s'%var

    histos = {}

    filename_nominal = 'quantiles/data/%s/plots/plotter.root'%(var)
    #filename_nominal = 'results/data/%s/plots/plotter.root'%(var)
    #filename_nominal = 'results/data/%s/plots/%s_unfolded.root'%(var,var)
    
    fIn_n=ROOT.TFile.Open(filename_nominal)

    recList_nominal  = getPathToObjects(fIn_n.Get('%s_rec'%(var)))
    genList_nominal  = getPathToObjects(fIn_n.Get('%s_gen'%(var)))

    if compare == 'gen':
      for p in genList_nominal:
          h=fIn_n.Get(p)
          hname=h.GetName()
          if 'TTJets' in hname or 'SingleT' in hname:
            if not '172' in histos : 
                histos['172']=h.Clone('172')
            else : 
               histos['172'].Add(h)
    
    if compare == 'rec':
      for p in recList_nominal:
          h=fIn_n.Get(p)
          hname=h.GetName()
          if 'TTJets' in hname or 'SingleT' in hname:
            if not '172' in histos : 
                histos['172']=h.Clone('172')
            else : 
               histos['172'].Add(h)

    for m in masspoints:

        # Open plotter.root for each of the mass files which contains all of the histograms
        #filename = 'results/mc/%s/%d/plots/plotter.root'%(var,m)
        filename = 'quantiles/mc/%s/%d/plots/plotter.root'%(var,m)

        m = str(m)

        fIn=ROOT.TFile.Open(filename)

        recList       = getPathToObjects(fIn.Get('%s_rec'%(var)))
        genList       = getPathToObjects(fIn.Get('%s_gen'%(var)))

        if compare == 'gen':
          for p in genList:
            h=fIn.Get(p)
            hname=h.GetName()
            if 'TTJets' in hname:
              if not m in histos : 
                  histos[m]=h.Clone(m)
              else : 
                 histos[m].Add(h)

        if compare == 'rec':
          for p in recList:
            h=fIn.Get(p)
            hname=h.GetName()
            if 'TTJets' in hname:
              if not m in histos : 
                  histos[m]=h.Clone(m)
              else : 
                 histos[m].Add(h)

    ratplot = RatioPlot('ratioplot')
    ratplot.ratiotitle = "Ratio wrt 172.5 GeV"
    ratplot.ratiorange = (0.5, 1.5)

    reference = [histos['172']]
    ratplot.reference = reference

    colors = [ROOT.kViolet-6,
              ROOT.kBlue-3, 
              ROOT.kRed-4, 
              ROOT.kOrange-3,
              ROOT.kGreen-3,
              ROOT.kGreen+2,
              ROOT.kAzure-2,]
    
    ratplot.colors1=colors
    ratplot.legpos = (0.68, 0.40)
    #ratplot.legpos = (0.64, 0.10)
    ratplot.tagpos = (0.88,0.87)
    #ratplot.normalized = False

    for h in sorted(histos.keys()):
        legentry = 'm_{t} = %s.5 GeV' % h
        histos[h].SetFillColor(0)
        try:
          histo = histos[h]
          ratplot.add(histo, legentry)
        except KeyError: pass
   
    if compare == 'gen':
      ratplot.tag = 'Generated level'
      ratplot.show("%s_gen_comparisonNorm"%var, outDir)
    if compare == 'rec':
      ratplot.tag = 'Reconstructed level'
      ratplot.show("%s_rec_comparisonNorm"%var, outDir)
    
def graphMass(v):

     var = str(v)
     #filename = 'results/mc/'+var+'/moments_'+var+'_unf.txt'
     #filename_gen = 'results/mc/'+var+'/moments_'+var+'_gen.txt'
     filename = 'quantiles/mc/'+var+'/moments_'+var+'_unf.txt'
     filename_gen = 'quantiles/mc/'+var+'/moments_'+var+'_gen.txt'
     filename_rec = 'quantiles/mc/'+var+'/moments_'+var+'_rec.txt'

     #pT positive lepton
     if var == 'ptpos': title = 'p_{T}(l^{+})'
  
     #pT charged-lepton pair
     if var == 'ptll': title = 'p_{T}(l^{+}l^{-})'

     #M charged-lepton pair
     if var == 'mll': title = 'M(l^{+}l^{-})'
        
     #Scalar sum of E
     if var == 'EposEm': title = 'E(l^{+})+E(l^{-})'
       
     #Scalar sum of Pt
     if var == 'ptposptm': title = 'p_{T}(l^{+})+p_{T}(l^{-})'

     m_unf,u0_unf,u1_unf,erru1_unf,u2_unf,erru2_unf = numpy.loadtxt(filename,unpack =True)
     m_gen,u0_gen,u1_gen,erru1_gen,u2_gen,erru2_gen = numpy.loadtxt(filename_gen,unpack =True)
     m_rec,u0_rec,u1_rec,erru1_rec,u2_rec,erru2_rec = numpy.loadtxt(filename_rec,unpack =True)

     mass_rec = array('d', m_rec)
     u0_rec = array('d', u0_rec)
     e0_rec = array('d', [0.0]*len(mass_rec))
     u1_rec = array('d', u1_rec)
     erru1_rec = array('d', erru1_rec)
     u2_rec = array('d', u2_rec)
     erru2_rec = array('d', erru2_rec)

     mass_gen= array('d', m_gen)
     u0_gen = array('d', u0_gen)
     e0_gen = array('d', [0.0]*len(mass_gen))
     u1_gen = array('d', u1_gen)
     erru1_gen = array('d', erru1_gen)
     u2_gen = array('d', u2_gen)
     erru2_gen = array('d', erru2_gen)

     mass_unf= array('d', m_unf)
     u0_unf = array('d', u0_unf)
     e0_unf = array('d', [0.0]*len(mass_unf))
     u1_unf = array('d', u1_unf)
     erru1_unf = array('d', erru1_unf)
     u2_unf = array('d', u2_unf)
     erru2_unf = array('d', erru2_unf)

     graphs = {}

     graphs['u1_'+var+'_unf'] = ROOT.TGraphErrors(len(mass_unf),mass_unf,u1_unf,e0_unf,erru1_unf)
     graphs['u2_'+var+'_unf'] = ROOT.TGraphErrors(len(mass_unf),mass_unf,u2_unf,e0_unf,erru2_unf)

     graphs['u1_'+var+'_rec'] = ROOT.TGraphErrors(len(mass_rec),mass_rec,u1_rec,e0_rec,erru1_rec)
     graphs['u2_'+var+'_rec'] = ROOT.TGraphErrors(len(mass_rec),mass_rec,u2_rec,e0_rec,erru2_rec)

     graphs['u1_'+var+'_gen'] = ROOT.TGraphErrors(len(mass_gen),mass_gen,u1_gen,e0_gen,erru1_gen)
     graphs['u2_'+var+'_gen'] = ROOT.TGraphErrors(len(mass_gen),mass_gen,u2_gen,e0_gen,erru2_gen)

     graphs['u1_'+var] = ROOT.TGraphErrors(len(u1_unf),u1_unf,u1_gen,erru1_unf,erru1_gen)
     graphs['u2_'+var] = ROOT.TGraphErrors(len(u2_unf),u2_unf,u2_gen,erru2_unf,erru2_gen)

     graphs['u1_'+var+'_r'] = ROOT.TGraphErrors(len(u1_rec),u1_rec,u1_gen,erru1_rec,erru1_gen)
     graphs['u2_'+var+'_r'] = ROOT.TGraphErrors(len(u2_rec),u2_rec,u2_gen,erru2_rec,erru2_gen)

     for g in graphs:
          gStyle.SetOptStat(0)
          gStyle.SetOptTitle(0)

          outputName = 'quantiles/mc/'+var+'/'+g
          
          #Create a canvas for plotting your graph
          c = ROOT.TCanvas(g,g)
          gStyle.SetCanvasDefH(600);
          gStyle.SetCanvasDefW(800);
          c.SetLeftMargin(0.15);
          c.SetRightMargin(0.25)
          c.SetBottomMargin(0.25);
          p1=ROOT.TPad('p1','p1',0.,0.,1.0,1.0)
          p1.Draw()

          #Create CMS labels
          tlat1 = TLatex()
          tlat1.SetNDC()
          tlat1.SetTextFont(42)
          tlat1.SetTextSize(0.04)
          tlat1.SetTextAlign(31)
          tlat1.DrawLatex(0.16, 0.91,'#bf{CMS}')

          tlat2 = TLatex()
          tlat2.SetNDC()
          tlat2.SetTextFont(42)
          tlat2.SetTextSize(0.04)
          tlat2.SetTextAlign(31)
          tlat2.DrawLatex(0.30, 0.91, '#it{Simulation}')
          #tlat2.DrawLatex(0.38, 0.91, '#it{Simulation}')
          #tlat2.DrawLatex(0.81, 0.91, '19.7 fb^{-1}')
          #tlat2.DrawLatex(0.90, 0.91,'(8 TeV)')
          tlat2.DrawLatex(0.90, 0.91,'8 TeV')

          tleg = ROOT.TLegend(0.12, 0.75, .50, 0.89)
          tleg.SetBorderSize(0)
          tleg.SetFillColor(0)
          tleg.SetFillStyle(0)
          tleg.SetShadowColor(0)
          tleg.SetTextFont(43)
          tleg.SetTextSize(20)

          p1.cd()
          graphs[g].SetTitle("")
          graphs[g].SetMarkerColor(kBlack)
          graphs[g].SetMarkerStyle(8)
          graphs[g].SetMarkerSize(0.7)
          if g == 'u1_'+var+'_unf':
               graphs[g].GetXaxis().SetTitle("Top Mass [GeV]")
               graphs[g].GetYaxis().SetTitle("Unfolded #mu^{(1)}_{%s} [GeV]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(1)}_{%s}'%title, 'pl')
          if g == 'u1_'+var+'_gen':
               graphs[g].GetXaxis().SetTitle("Top Mass [GeV]")
               graphs[g].GetYaxis().SetTitle("Generated level #mu^{(1)}_{%s} [GeV]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(1)}_{%s}'%title, 'pl')
          if g == 'u1_'+var+'_rec':
               graphs[g].GetXaxis().SetTitle("Top Mass [GeV]")
               graphs[g].GetYaxis().SetTitle("Reconstructed level #mu^{(1)}_{%s} [GeV]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(1)}_{%s}'%title, 'pl')
          if g == 'u2_'+var+'_unf':
               graphs[g].GetXaxis().SetTitle("Top Mass [GeV]")
               graphs[g].GetYaxis().SetTitle("Unfolded #mu^{(2)}_{%s} [GeV^{2}]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(2)}_{%s}'%title, 'pl')
          if g == 'u2_'+var+'_gen':
               graphs[g].GetXaxis().SetTitle("Top Mass [GeV]")
               graphs[g].GetYaxis().SetTitle("Generated level #mu^{(2)}_{%s} [GeV^{2}]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(2)}_{%s}'%title, 'pl')
          if g == 'u2_'+var+'_rec':
               graphs[g].GetXaxis().SetTitle("Top Mass [GeV]")
               graphs[g].GetYaxis().SetTitle("Reconstructed level #mu^{(2)}_{%s} [GeV^{2}]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(2)}_{%s}'%title, 'pl')
          if g == 'u1_'+var:
               graphs[g].GetYaxis().SetTitle("Generated level #mu^{(1)}_{%s} [GeV]"%title)
               graphs[g].GetXaxis().SetTitle("Unfolded #mu^{(1)}_{%s} [GeV]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(1)}_{%s}'%title, 'pl')
          if g == 'u1_'+var+'_r':
               graphs[g].GetYaxis().SetTitle("Generated level #mu^{(1)}_{%s} [GeV]"%title)
               graphs[g].GetXaxis().SetTitle("Reconstructed #mu^{(1)}_{%s} [GeV]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(1)}_{%s}'%title, 'pl')
          if g == 'u2_'+var:
               graphs[g].GetYaxis().SetTitle("Generated level #mu^{(2)}_{%s} [GeV^{2}]"%title)
               graphs[g].GetXaxis().SetTitle("Unfolded #mu^{(2)}_{%s} [GeV^{2}]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(2)}_{%s}'%title, 'pl')
          if g == 'u2_'+var+ '_r':
               graphs[g].GetYaxis().SetTitle("Generated level #mu^{(2)}_{%s} [GeV^{2}]"%title)
               graphs[g].GetXaxis().SetTitle("Reconstructed #mu^{(2)}_{%s} [GeV^{2}]"%title)
               tleg.AddEntry(graphs[g], '#mu^{(2)}_{%s}'%title, 'pl')
          graphs[g].GetXaxis().SetTitleSize(0.037)
          graphs[g].GetXaxis().SetLabelSize(0.037)
          graphs[g].GetYaxis().SetTitleOffset(1.2)
          graphs[g].GetYaxis().SetTitleSize(0.037)
          graphs[g].GetYaxis().SetLabelSize(0.037)
          graphs[g].Draw("AP")

          #if g == 'u1_'+var or g == 'u2_'+var :
          #     f1 = ROOT.TF1("f1","x",0,5000);
          #     f1.SetLineColor(12)
          #     f1.SetLineWidth(1)
          #     f1.SetLineStyle(2)
          #     f1.Draw("same")
          #     tleg.AddEntry("f1", "Expected", "l")

          #Define the fit function and make the linear fit 
          tf = ROOT.TF1("linear", "pol1", 160., 190.)
          graphs[g].Fit(tf, "WQ")
          tf.SetLineColor(ROOT.kBlue)
          tf.SetLineWidth(2)
          tf.SetLineStyle(1)
          tf.Draw("same")
          
          tleg.Draw()

          #Get parameters from fit
          """
          p0=linear.GetParameter(0)
          p0Err=linear.GetParError(0)
          p1=linear.GetParameter(1)
          p1Err=linear.GetParError(1)
          chi2=linear.GetChisquare()
          NDF=linear.GetNDF()
          linear = graphs[g].GetFunction("linear")  
          """        

          #save and delete
          for ext in ['png','pdf']:
            c.SaveAs('%s.%s'%(outputName,ext))
          del c

def plotter(opt):

     var = ['ptpos']
     #var = ['ptpos','ptll','mll','EposEm','ptposptm']

     for v in var:

      #os.system("python scripts/runDileptonUnfolding.py -i /store/cmst3/group/top/summer2015/treedir_bbbcb36/ttbar/mass_scan -o quantiles/mc/%s/ --jobs 8 -m True -v %s"%(v,v))
           
      # Create an array with the different masses
      mass = [166,169,171,173,175,178]
      #mass = [166]
           
      opt.var = v
      
      for m in mass:
        m = str(m)

        # Get all the plots
        #os.system("python scripts/runPlotter.py -j test/topss2014/mass_scan_samples.json -o results/mc/%s/%s/plots results/mc/%s/%s"%(v,m,v,m))
        #os.system("python scripts/runDileptonUnfolding.py -r results/mc/%s/%s/plots/plotter.root -v %s -o results/mc/%s/%s"%(v,m,v,v,m)) 
        os.system("python scripts/runPlotter.py -j test/topss2014/mass_scan_samples.json -o quantiles/mc/%s/%s/plots quantiles/mc/%s/%s"%(v,m,v,m))
        os.system("python scripts/runDileptonUnfolding.py -r quantiles/mc/%s/%s/plots/plotter.root -v %s -o quantiles/mc/%s/%s"%(v,m,v,v,m)) 

        root = 'quantiles/mc/%s/%s/plots/%s_unfolded.root'%(v,m,v)
        getMoments(root,v,'unf',m)
        root = 'quantiles/mc/%s/%s/plots/plotter.root'%(v,m)
        getMoments(root,v,'rec',m)
        getMoments(root,v,'gen',m)

        #root = 'results/mc/%s/%s/plots/plotter.root'%(v,m)
        root = 'quantiles/mc/%s/%s/plots/plotter.root'%(v,m)

        getStability(root,v,m)

      graphMass(v)

      compare = 'gen'
      compareMass(compare,v)
      compare = 'rec'
      compareMass(compare,v)

def main():
     usage = 'usage: %prog [options]'
     parser = optparse.OptionParser(usage)
     parser.add_option('-g',
                          dest='gen',   
                          default='a', 
                          help='Create graphs for unfolded gen distributions depending on the mass [default: %default]')
     parser.add_option('-v', '--var',
                          dest='var', 
                          default='ptpos',
                          help='Variable to unfold (note requires var_rec,var_gen,var_migration plots stored in the ROOT file [default: %default]')
     parser.add_option('-o', '--output',
                          dest='output', 
                          default='unfoldResults',                                                                       
                          help='Output directory [default: %default]')
     parser.add_option('-r', '--root',
                          dest='root', 
                          default=None,
                          help='ROOT file with distributions to get moments[default: %default]')
     parser.add_option('-c',
                          dest='compare', 
                          default='gen',
                          help='Distributions to compare: gen, rec [default: %default]')
     parser.add_option('-m',
                          dest='moments', 
                          default='gen',
                          help='Distributions to get moments: gen, rec or unf [default: %default]')

     (opt, args) = parser.parse_args()

     ROOT.gStyle.SetOptStat(0)
     ROOT.gStyle.SetOptTitle(0)
     ROOT.gROOT.SetBatch(True)

     plotter(opt)

if __name__ == "__main__":
    sys.exit(main())
