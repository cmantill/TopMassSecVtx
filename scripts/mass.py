#!/usr/bin/env python

import ROOT
from ROOT import *
import os,sys,numpy
import array
import pprint

def main(argv=None):

     var = sys.argv[1]

     m,u0,u1,u2,erru1,erru2 = numpy.loadtxt('/afs/cern.ch/work/c/cmantill/public/CMSSW_5_3_22/src/UserCode/TopMassSecVtx/results/'+var+'/moments_'+var+'.txt',unpack =True)
     #m_gen,u0_gen,u1_gen,u2_gen,erru1_gen,erru2_gen = numpy.loadtxt('/afs/cern.ch/work/c/cmantill/public/CMSSW_5_3_22/src/UserCode/TopMassSecVtx/results/'+var+'/momentsgen_'+var+'.txt',unpack =True)

     mass = array.array('d', m)
     u0 = array.array('d', u0)
     e0 = array.array('d', [0.0]*len(mass))
     u1 = array.array('d', u1)
     erru1 = array.array('d', erru1)
     u2 = array.array('d', u2)
     erru2 = array.array('d', erru2)

     graphs = {}

     graphs['u0_'+var] = ROOT.TGraph(len(mass),mass,u0)
     graphs['u1_'+var] = ROOT.TGraphErrors(len(mass),mass,u1,e0,erru1)
     graphs['u2_'+var] = ROOT.TGraphErrors(len(mass),mass,u2,e0,erru2)

     for g in graphs:
          gStyle.SetOptStat(0)
          gStyle.SetOptTitle(0)

          #Create a canvas for plotting your graph
          c = ROOT.TCanvas('c','c')
          gStyle.SetCanvasDefH(400);
          gStyle.SetCanvasDefW(700);
          gStyle.SetPadLeftMargin(0.01);
          gStyle.SetPadRightMargin(0.01);    

          #Define the fit function and make the linear fit 
          linear = ROOT.TF1("linear", "pol1",160,180)
          linear.SetLineColor(kBlue)
          linear.SetLineWidth(2)
          linear.SetLineStyle(1) 
          graphs[g].Fit('pol1')

          graphs[g].SetTitle("")
          graphs[g].SetLineColor(2);
          graphs[g].SetLineWidth(2)
          graphs[g].SetMarkerColor(4)
          graphs[g].SetMarkerStyle(8)
          graphs[g].GetXaxis().SetTitle("Top Mass [GeV/c^{2}]")
          graphs[g].GetXaxis().SetTitleSize(0.047)
          graphs[g].GetXaxis().SetLabelSize(0.047)
          if g == 'u0_'+var:
               graphs[g].GetYaxis().SetTitle("Mellin Moment u_{0}=1")
          if g == 'u1_'+var:
               graphs[g].GetYaxis().SetTitle("First Mellin Moment u_{1} = #LT p_{T}(l^{+}) #GT")
          if g == 'u2_'+var:
               graphs[g].GetYaxis().SetTitle("Second Mellin Moment u_{2} = #LT p_{T}(l^{+}) #GT^{2}")
          graphs[g].GetYaxis().SetTitleSize(0.047)
          graphs[g].GetYaxis().SetLabelSize(0.047)
          graphs[g].Draw()

          #Get parameters from fit
          p0=linear.GetParameter(0)
          p0Err=linear.GetParError(0)
          p1=linear.GetParameter(1)
          p1Err=linear.GetParError(1)
          chi2=linear.GetChisquare()
          NDF=linear.GetNDF()
          linear = graphs[g].GetFunction("linear")
     
          #Create some labels about the statistics
          #caption=ROOT.TLatex()
          #caption.SetTextSize(0.037)
          #caption.SetTextFont(30)
          #caption.SetNDC()
          #caption.DrawLatex(0.1,0.94,'#bf{CMS Work in Progress, #sqrt{s}=8TeV}')
          #caption.DrawLatex(0.75,0.84,'Fit Results')
          #caption.DrawLatex(0.7,0.79,'y=%3.4f + %3.4f x #pm %3.4f %3.4f'%(p0,p1,p0Err,p1Err))
          #caption.DrawLatex(0.7,0.69,'#chi^{2}/ndf=%3.4f/%3.4f'%(chi2,NDF))

          #save and delete
          outputName = '/afs/cern.ch/work/c/cmantill/public/CMSSW_5_3_22/src/UserCode/TopMassSecVtx/results/'+var+'/'+g+'_'+var
          c.SaveAs(outputName+'.pdf')
          del c
          #del caption

     filename = '/afs/cern.ch/work/c/cmantill/public/CMSSW_5_3_22/src/UserCode/TopMassSecVtx/results/'+var+'/moments_'+var+'.root'
     fOut=ROOT.TFile.Open(filename,'RECREATE')
     for g in graphs: graphs[g].Write(g)
     print 'Graphs saved in %s' % fOut.GetName()
     fOut.Close()

if __name__ == "__main__":
    sys.exit(main())
