#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TList.h"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include <stdio.h>
#include <TROOT.h>
#include "TSystem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>


void savePlot(TCanvas*, TString);

void RedPercent_2D_StaMuonHits_2015() {

    // Title
    // -----
    std::string title = "62Xvs71X";


    // Main Comparison
    // ---------------
    std::string label1 = "CMSSW 62X SLHC (2023 Geometry)"; 
    std::string label2 = "CMSSW 71X (2015 Geometry)"; 

    // Additional Information
    // ----------------------
    std::string info1 = "Hits used in Stand Alone Muon Fit";
    std::string info2 = "these are Segments (CSC,DT) and Rechits (RPC)";
    std::string info3 = "Sample used: RelValZmumuJets_Pt_20_300";
    // std::string info4 = "Muon Selection :: only STA Muons matched to Z boson";
    std::string info4 = "Muon Selection :: Tight Muon ID + p_{T} > 10 GeV/c";
    // std::string info4 = "Muon Selection :: Loose Muon ID + p_{T} > 10 GeV/c";

    double y_range_1 = 00.00;
    double y_range_2 = 15.00;
    double x_title_offset = 0.75;
    double y_title_offset = 0.75;
    double z_title_offset = 0.75;

    std::string xaxislabel = "#eta_{STA Muon}";
    std::string yaxislabel = "\# Hits used in fit";
    std::string zaxislabel = "Percentage";

    // Reading the filename from command line
    // --------------------------------------
    // std::string RootFileName = gFile->GetName();
    
    // Filename hardcoded
    // ------------------
    std::string RootFileName1 = "STAMuon_62X_2023.root";
    std::string RootFileName2 = "STAMuon_70X_2015.root";

    // Filename base
    // -------------
    std::string Histo_All   = "SegmentsAndHits/Rechits_All_Eta_2D_Reduced";
    std::string Histo_RPC   = "SegmentsAndHits/Rechits_RPC_Eta_2D_Reduced";
    std::string Histo_CSC   = "SegmentsAndHits/Rechits_CSC_Eta_2D_Reduced";
    std::string Histo_DT    = "SegmentsAndHits/Rechits_DT_Eta_2D_Reduced";
    std::string Histo_Denom = "Muon_All_Reduced";

    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1,0);
    gStyle->SetOptStat(000000);
    gStyle->SetOptTitle(0);
    gStyle->SetGridStyle(1);    // grid
    gROOT->UseCurrentStyle();
    
    TFile * infile1             = new TFile((RootFileName1).c_str(),"READ");     std::cout<<"File 1  = "<<infile1<<std::endl; 
    TFile * infile2             = new TFile((RootFileName2).c_str(),"READ");     std::cout<<"File 2  = "<<infile2<<std::endl; 

    TH2F * All_histo1           = (TH2F*) infile1->Get(Histo_All.c_str());   std::cout<<"Histo  = "<<All_histo1<<std::endl; 
    TH2F * All_histo2           = (TH2F*) infile2->Get(Histo_All.c_str());   std::cout<<"Histo  = "<<All_histo2<<std::endl; 
    TH2F * RPC_histo1           = (TH2F*) infile1->Get(Histo_RPC.c_str());   std::cout<<"Histo  = "<<RPC_histo1<<std::endl; 
    TH2F * RPC_histo2           = (TH2F*) infile2->Get(Histo_RPC.c_str());   std::cout<<"Histo  = "<<RPC_histo2<<std::endl; 
    TH2F * CSC_histo1           = (TH2F*) infile1->Get(Histo_CSC.c_str());   std::cout<<"Histo  = "<<CSC_histo1<<std::endl; 
    TH2F * CSC_histo2           = (TH2F*) infile2->Get(Histo_CSC.c_str());   std::cout<<"Histo  = "<<CSC_histo2<<std::endl; 
    TH2F * DT_histo1            = (TH2F*) infile1->Get(Histo_DT.c_str());    std::cout<<"Histo  = "<<DT_histo1<<std::endl; 
    TH2F * DT_histo2            = (TH2F*) infile2->Get(Histo_DT.c_str());    std::cout<<"Histo  = "<<DT_histo2<<std::endl; 

    TH1F * Denom_histo1         = (TH1F*) infile1->Get(Histo_Denom.c_str());    std::cout<<"Histo  = "<<Denom_histo1<<std::endl; 
    TH1F * Denom_histo2         = (TH1F*) infile2->Get(Histo_Denom.c_str());    std::cout<<"Histo  = "<<Denom_histo2<<std::endl; 

    All_histo1->SetLineColor(kGray+2); All_histo1->SetLineWidth(1);         All_histo2->SetLineColor(kRed); All_histo2->SetLineWidth(1); 
    All_histo1->SetFillColor(kBlue-7); /*All_histo1->SetFillStyle(4050);*/  All_histo2->SetFillColor(kRed); All_histo2->SetFillStyle(3002); 
    RPC_histo1->SetLineColor(kGray+2); RPC_histo1->SetLineWidth(1);         RPC_histo2->SetLineColor(kRed); RPC_histo2->SetLineWidth(1); 
    RPC_histo1->SetFillColor(kBlue-7); /*RPC_histo1->SetFillStyle(4050);*/  RPC_histo2->SetFillColor(kRed); RPC_histo2->SetFillStyle(3002); 
    CSC_histo1->SetLineColor(kGray+2);  CSC_histo1->SetLineWidth(1);          CSC_histo2->SetLineColor(kRed);  CSC_histo2->SetLineWidth(1); 
    CSC_histo1->SetFillColor(kBlue-7);  /*CSC_histo1->SetFillStyle(4050);*/   CSC_histo2->SetFillColor(kRed);  CSC_histo2->SetFillStyle(3002); 
    DT_histo1->SetLineColor(kGray+2);   DT_histo1->SetLineWidth(1);           DT_histo2->SetLineColor(kRed);   DT_histo2->SetLineWidth(1); 
    DT_histo1->SetFillColor(kBlue-7);   /*DT_histo1->SetFillStyle(4050);*/    DT_histo2->SetFillColor(kRed);   DT_histo2->SetFillStyle(3002); 
    // std::string Xaxis = "log_{10} TOF (ns)";
    std::string Xaxis = "#eta";
    std::string Yaxis = "\# Hits";

    All_histo1->SetTitle("All");  All_histo1->GetXaxis()->SetTitle(Xaxis.c_str()); All_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    RPC_histo1->SetTitle("RPC");  RPC_histo1->GetXaxis()->SetTitle(Xaxis.c_str()); RPC_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    CSC_histo1->SetTitle("CSC");   CSC_histo1->GetXaxis()->SetTitle(Xaxis.c_str());  CSC_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    DT_histo1->SetTitle("DT");     DT_histo1->GetXaxis()->SetTitle(Xaxis.c_str());   DT_histo1->GetYaxis()->SetTitle(Yaxis.c_str());

    double denominator, numerator;
    for(int i=0; i<Denom_histo2->GetXaxis()->GetNbins(); ++i) {
      denominator = Denom_histo2->GetBinContent(i+1);
      for(int j=0; j<All_histo2->GetYaxis()->GetNbins(); ++j) {
	numerator = All_histo2->GetBinContent(i+1,j+1);
	if (numerator > 0.0 && denominator > 0.0) All_histo2->SetBinContent(i+1,j+1, numerator/denominator*100); All_histo2->SetBinError(i+1,j+1, sqrt(numerator)/denominator*100);
	numerator = RPC_histo2->GetBinContent(i+1,j+1);
	if (numerator > 0.0 && denominator > 0.0) RPC_histo2->SetBinContent(i+1,j+1, numerator/denominator*100); RPC_histo2->SetBinError(i+1,j+1, sqrt(numerator)/denominator*100);
	numerator = CSC_histo2->GetBinContent(i+1,j+1);
	if (numerator > 0.0 && denominator > 0.0) CSC_histo2->SetBinContent(i+1,j+1, numerator/denominator*100);  CSC_histo2->SetBinError(i+1,j+1, sqrt(numerator)/denominator*100);
	numerator = DT_histo2->GetBinContent(i+1,j+1);
	if (numerator > 0.0 && denominator > 0.0) DT_histo2->SetBinContent(i+1,j+1, numerator/denominator*100);   DT_histo2->SetBinError(i+1,j+1, sqrt(numerator)/denominator*100);
      }
    }


    int All_1_entries = All_histo1->GetEntries();       int All_2_entries = All_histo2->GetEntries();
    int RPC_1_entries = RPC_histo1->GetEntries();       int RPC_2_entries = RPC_histo2->GetEntries();
    int CSC_1_entries  = CSC_histo1->GetEntries();        int CSC_2_entries  = CSC_histo2->GetEntries();
    int DT_1_entries   = DT_histo1->GetEntries();         int DT_2_entries   = DT_histo2->GetEntries();

    std::stringstream ss_All_1_label; ss_All_1_label << label1;      std::stringstream ss_All_2_label; ss_All_2_label << label2;
    std::stringstream ss_RPC_1_label; ss_RPC_1_label << label1;      std::stringstream ss_RPC_2_label; ss_RPC_2_label << label2;
    std::stringstream ss_CSC_1_label;   ss_CSC_1_label   << label1;      std::stringstream ss_CSC_2_label;   ss_CSC_2_label   << label2;
    std::stringstream ss_DT_1_label;    ss_DT_1_label    << label1;      std::stringstream ss_DT_2_label;    ss_DT_2_label    << label2;

    std::string All_1_label = ss_All_1_label.str();    std::string All_2_label = ss_All_2_label.str();
    std::string RPC_1_label = ss_RPC_1_label.str();    std::string RPC_2_label = ss_RPC_2_label.str();
    std::string CSC_1_label   = ss_CSC_1_label.str();      std::string CSC_2_label = ss_CSC_2_label.str();
    std::string DT_1_label    = ss_DT_1_label.str();       std::string DT_2_label = ss_DT_2_label.str();

    /*
    All_histo1->Scale(1.0/All_1_entries); All_histo2->Scale(1.0/All_2_entries);
    RPC_histo1->Scale(1.0/RPC_1_entries); RPC_histo2->Scale(1.0/RPC_2_entries);
    CSC_histo1->Scale(1.0/CSC_1_entries);   CSC_histo2->Scale(1.0/CSC_2_entries);
    DT_histo1->Scale(1.0/DT_1_entries);     DT_histo2->Scale(1.0/DT_2_entries);
    */

    double leg_x1 = 0.55, leg_x2 = 0.85, leg_y1 = 0.700, leg_y2 = 0.850; 

    TLegend * l1 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l1->SetLineColor(1);    l1->SetLineStyle(1);  l1->SetLineWidth(1);  l1->SetFillColor(4000); l1->SetBorderSize(1);
    TLegend * l2 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l2->SetLineColor(1);    l2->SetLineStyle(1);  l2->SetLineWidth(1);  l2->SetFillColor(4000); l2->SetBorderSize(1);
    TLegend * l3 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l3->SetLineColor(1);    l3->SetLineStyle(1);  l3->SetLineWidth(1);  l3->SetFillColor(4000); l3->SetBorderSize(1);
    TLegend * l4 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l4->SetLineColor(1);    l4->SetLineStyle(1);  l4->SetLineWidth(1);  l4->SetFillColor(4000); l4->SetBorderSize(1);

    l1->AddEntry(All_histo1, All_1_label.c_str(), "F");     l1->AddEntry(All_histo2, All_2_label.c_str(), "F");
    l2->AddEntry(RPC_histo1, RPC_1_label.c_str(), "F");     l2->AddEntry(RPC_histo2, RPC_2_label.c_str(), "F");
    l3->AddEntry(CSC_histo1,  CSC_1_label.c_str(),   "F");     l3->AddEntry(CSC_histo2,  CSC_2_label.c_str(),   "F");
    l4->AddEntry(DT_histo1,   DT_1_label.c_str(),    "F");     l4->AddEntry(DT_histo2,   DT_2_label.c_str(),    "F");

    TLatex latex_cmslab; latex_cmslab.SetNDC(); latex_cmslab.SetTextSize(0.035); latex_cmslab.SetTextAlign(11);
    TLatex latex_info;   latex_info.SetNDC();   latex_info.SetTextSize(0.025);    latex_info.SetTextAlign(31);
    TLatex latex_infoBig;latex_infoBig.SetNDC();latex_infoBig.SetTextSize(0.035); latex_infoBig.SetTextAlign(31);

    double infoX_x_right = 0.6275; 
    double info1_y = 0.825, info2_y = 0.775, info3_y = 0.725, info4_y = 0.675;

    /*
    TLegend * l5 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l5->SetLineColor(1);    l5->SetLineStyle(1);  l5->SetLineWidth(1);  l5->SetFillColor(4000); l5->SetBorderSize(1);
    TLegend * l6 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l6->SetLineColor(1);    l6->SetLineStyle(1);  l6->SetLineWidth(1);  l6->SetFillColor(4000); l6->SetBorderSize(1);
    TLegend * l7 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l7->SetLineColor(1);    l7->SetLineStyle(1);  l7->SetLineWidth(1);  l7->SetFillColor(4000); l7->SetBorderSize(1);
    TLegend * l8 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l8->SetLineColor(1);    l8->SetLineStyle(1);  l8->SetLineWidth(1);  l8->SetFillColor(4000); l8->SetBorderSize(1);

    l5->AddEntry(All_histo1,  All_1_label.c_str(), "F");       l5->AddEntry(All_histo2,  All_2_label.c_str(), "F");
    l6->AddEntry(RPC_histo1,  RPC_1_label.c_str(), "F");       l6->AddEntry(RPC_histo2,  RPC_2_label.c_str(), "F");
    l7->AddEntry(CSC_histo1,  CSC_1_label.c_str(),   "F");     l7->AddEntry(CSC_histo2,  CSC_2_label.c_str(),   "F");
    l8->AddEntry(DT_histo1,   DT_1_label.c_str(),    "F");     l8->AddEntry(DT_histo2,   DT_2_label.c_str(),    "F");
    */

    // Save the Other Canvasses
    y_range_2 = 15.00;
    std::stringstream pdfname2ss; pdfname2ss<<"STAMuonHits_2D_2015_"<<"RedPercent"<<"_All"; std::string pdfname2 = pdfname2ss.str(); 
    TCanvas * c2 = new TCanvas(pdfname2.c_str(), "", 800, 400); c2->cd(); c2->SetTicks(1,1); c2->SetLeftMargin(0.05); c2->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    All_histo2->SetMarkerSize(1.75); All_histo2->Draw("textcolz"); All_histo2->GetYaxis()->SetRangeUser(y_range_1,y_range_2); All_histo2->GetZaxis()->SetRangeUser(0,100); 
    All_histo2->GetXaxis()->SetTitleOffset(x_title_offset); All_histo2->GetYaxis()->SetTitleOffset(y_title_offset);   All_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    All_histo2->GetXaxis()->SetTitle(xaxislabel.c_str());   All_histo2->GetYaxis()->SetTitle(yaxislabel.c_str());     All_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "All");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); 
    latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str()); latex_info.DrawLatex(infoX_x_right, info4_y, info4.c_str());
    savePlot(c2, pdfname2.c_str());

    y_range_2 = 10.00;
    std::stringstream pdfname3ss; pdfname3ss<<"STAMuonHits_2D_2015_"<<"RedPercent"<<"_RPC"; std::string pdfname3 = pdfname3ss.str(); 
    TCanvas * c3 = new TCanvas(pdfname3.c_str(), "", 800, 400); c3->cd(); c3->SetTicks(1,1); c3->SetLeftMargin(0.05); c3->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    RPC_histo2->SetMarkerSize(1.75); RPC_histo2->Draw("textcolz"); RPC_histo2->GetYaxis()->SetRangeUser(y_range_1,y_range_2); RPC_histo2->GetZaxis()->SetRangeUser(0,100);
    RPC_histo2->GetXaxis()->SetTitleOffset(x_title_offset); RPC_histo2->GetYaxis()->SetTitleOffset(y_title_offset);   RPC_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    RPC_histo2->GetXaxis()->SetTitle(xaxislabel.c_str());   RPC_histo2->GetYaxis()->SetTitle(yaxislabel.c_str());     RPC_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "RPC");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c3, pdfname3.c_str());

    y_range_2 = 10.00;
    std::stringstream pdfname4ss; pdfname4ss<<"STAMuonHits_2D_2015_"<<"RedPercent"<<"_DT"; std::string pdfname4 = pdfname4ss.str(); 
    TCanvas * c4 = new TCanvas(pdfname4.c_str(), "", 800, 400); c4->cd(); c4->SetTicks(1,1); c4->SetLeftMargin(0.05); c4->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    DT_histo2->SetMarkerSize(1.75); DT_histo2->Draw("textcolz"); DT_histo2->GetYaxis()->SetRangeUser(y_range_1,y_range_2); DT_histo2->GetZaxis()->SetRangeUser(0,100); 
    DT_histo2->GetXaxis()->SetTitleOffset(x_title_offset); DT_histo2->GetYaxis()->SetTitleOffset(y_title_offset);   DT_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    DT_histo2->GetXaxis()->SetTitle(xaxislabel.c_str());   DT_histo2->GetYaxis()->SetTitle(yaxislabel.c_str());     DT_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "DT");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c4, pdfname4.c_str());

    y_range_2 = 10.00;
    std::stringstream pdfname5ss; pdfname5ss<<"STAMuonHits_2D_2015_"<<"RedPercent"<<"_CSC"; std::string pdfname5 = pdfname5ss.str(); 
    TCanvas * c5 = new TCanvas(pdfname5.c_str(), "", 800, 400); c5->cd(); c5->SetTicks(1,1); c5->SetLeftMargin(0.05); c5->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    CSC_histo2->SetMarkerSize(1.75); CSC_histo2->Draw("textcolz"); CSC_histo2->GetYaxis()->SetRangeUser(y_range_1,y_range_2); CSC_histo2->GetZaxis()->SetRangeUser(0,100);
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "CSC");
    CSC_histo2->GetXaxis()->SetTitleOffset(x_title_offset); CSC_histo2->GetYaxis()->SetTitleOffset(y_title_offset);   CSC_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    CSC_histo2->GetXaxis()->SetTitle(xaxislabel.c_str());   CSC_histo2->GetYaxis()->SetTitle(yaxislabel.c_str());     CSC_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c5, pdfname5.c_str());

    // reduced eta binning
    double reduced[] = {-2.5, -1.7, -1.6, -1.1, -0.9, -0.8, -0.3, -0.2, 0.2, 0.3, 0.8, 0.9, 1.1, 1.6, 1.7, 2.5};
    TH1F * Ineff_All_histo = new TH1F("Ineff_All_histo", "Ineff_All_histo", 15, reduced);
    for(int i = 0; i<15; ++i) {
      double ineff_value = All_histo2->GetBinContent(i+1,2+1);  double ineff_uncrt = All_histo2->GetBinError(i+1,2+1);
      std::cout<<"bin ["<<i+1<<"] inefficiency = "<<ineff_value<<" +/- "<<ineff_uncrt<<std::endl;
      Ineff_All_histo->SetBinContent(i+1, 100-(ineff_value/1.00)); Ineff_All_histo->SetBinError(i+1, ineff_uncrt/1.00);
    }

    y_range_1 = 90.00;
    y_range_2 = 105.00;
    std::string yaxislabel_ineff = "Efficiency";
    std::string info5 = "Stand Alone Muon Reco Efficiency";
    std::string info6 = "one RPC hit / DT,CSC segment missing";
    y_title_offset = 1.25;
    Ineff_All_histo->SetLineWidth(2.0); Ineff_All_histo->SetMarkerStyle(20); Ineff_All_histo->SetMarkerSize(1.25);
    std::stringstream pdfname6ss; pdfname6ss<<"STAMuonHits_2D_2015_"<<"RedPercent"<<"_Eff_1SegRemoved_All"; std::string pdfname6 = pdfname6ss.str();
    TCanvas * c6 = new TCanvas(pdfname6.c_str(), "", 800, 400); c6->cd(); c6->SetTicks(1,1); // c6->SetLeftMargin(0.05); c6->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    Ineff_All_histo->SetMarkerSize(1.75); Ineff_All_histo->Draw("PE"); Ineff_All_histo->GetYaxis()->SetRangeUser(y_range_1,y_range_2); // Ineff_All_histo->GetZaxis()->SetRangeUser(0,100);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "All");
    Ineff_All_histo->GetXaxis()->SetTitleOffset(x_title_offset); Ineff_All_histo->GetYaxis()->SetTitleOffset(y_title_offset);         // Ineff_All_histo->GetZaxis()->SetTitleOffset(z_title_offset);
    Ineff_All_histo->GetXaxis()->SetTitle(xaxislabel.c_str());   Ineff_All_histo->GetYaxis()->SetTitle(yaxislabel_ineff.c_str());     // Ineff_All_histo->GetZaxis()->SetTitle(zaxislabel.c_str());  
    latex_infoBig.DrawLatex(infoX_x_right, info1_y, info5.c_str()); latex_infoBig.DrawLatex(infoX_x_right, info2_y, info6.c_str()); 
    latex_infoBig.DrawLatex(infoX_x_right, info3_y, info3.c_str()); latex_infoBig.DrawLatex(infoX_x_right, info4_y, info4.c_str());
    savePlot(c6, pdfname6.c_str());



}

void savePlot(TCanvas * c, TString name) {
  c->SaveAs(name+".png");
  c->SaveAs(name+".eps");
  gSystem->Exec("epstopdf "+name+".eps");
  std::cout<<"Plot "<<name<<" saved"<<std::endl;
}
