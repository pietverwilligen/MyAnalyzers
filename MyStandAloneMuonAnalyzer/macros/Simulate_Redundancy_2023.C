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

int factorial(int);
void savePlot(TCanvas*, TString);

void Simulate_Redundancy_2023() {

    // Title
    // -----
    std::string title = "62Xvs71X";


    // Main Comparison
    // ---------------
    std::string label2 = "CMSSW 62X SLHC (2023 Geometry)"; 
    std::string label1 = "CMSSW 71X (2015 Geometry)"; 

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
    double x_title_offset = 1.00;
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
    std::string RootFileName2 = "STAMuon_62X_2023.root";
    std::string RootFileName1 = "STAMuon_70X_2015.root";

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
    // gStyle->SetPaintTextFormat("3.3f");
    // gStyle->SetTextColor(kWhite);
    std::stringstream pdfname2ss; pdfname2ss<<"STAMuonHits_2D_2023_"<<"RedPercent"<<"_All"; std::string pdfname2 = pdfname2ss.str(); 
    TCanvas * c2 = new TCanvas(pdfname2.c_str(), "", 800, 400); c2->cd(); c2->SetTicks(1,1); c2->SetLeftMargin(0.05); c2->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    All_histo2->SetMarkerSize(1.75); All_histo2->Draw("textcolz"); All_histo2->GetYaxis()->SetRangeUser(y_range_1,y_range_2); All_histo2->GetZaxis()->SetRangeUser(0,100); 
    All_histo2->GetXaxis()->SetTitleOffset(x_title_offset); All_histo2->GetYaxis()->SetTitleOffset(y_title_offset);   All_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    All_histo2->GetXaxis()->SetTitle(xaxislabel.c_str());   All_histo2->GetYaxis()->SetTitle(yaxislabel.c_str());     All_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "All");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); 
    latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str()); latex_info.DrawLatex(infoX_x_right, info4_y, info4.c_str());
    savePlot(c2, pdfname2.c_str());

    /*
    y_range_2 = 10.00;
    std::stringstream pdfname3ss; pdfname3ss<<"STAMuonHits_2D_2023_"<<"RedPercent"<<"_RPC"; std::string pdfname3 = pdfname3ss.str(); 
    TCanvas * c3 = new TCanvas(pdfname3.c_str(), "", 800, 400); c3->cd(); c3->SetTicks(1,1); c3->SetLeftMargin(0.05); c3->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    RPC_histo2->SetMarkerSize(1.75); RPC_histo2->Draw("textcolz"); RPC_histo2->GetYaxis()->SetRangeUser(y_range_1,y_range_2); RPC_histo2->GetZaxis()->SetRangeUser(0,100);
    RPC_histo2->GetXaxis()->SetTitleOffset(x_title_offset); RPC_histo2->GetYaxis()->SetTitleOffset(y_title_offset);   RPC_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    RPC_histo2->GetXaxis()->SetTitle(xaxislabel.c_str());   RPC_histo2->GetYaxis()->SetTitle(yaxislabel.c_str());     RPC_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "RPC");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c3, pdfname3.c_str());

    y_range_2 = 10.00;
    std::stringstream pdfname4ss; pdfname4ss<<"STAMuonHits_2D_2023_"<<"RedPercent"<<"_DT"; std::string pdfname4 = pdfname4ss.str(); 
    TCanvas * c4 = new TCanvas(pdfname4.c_str(), "", 800, 400); c4->cd(); c4->SetTicks(1,1); c4->SetLeftMargin(0.05); c4->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    DT_histo2->SetMarkerSize(1.75); DT_histo2->Draw("textcolz"); DT_histo2->GetYaxis()->SetRangeUser(y_range_1,y_range_2); DT_histo2->GetZaxis()->SetRangeUser(0,100); 
    DT_histo2->GetXaxis()->SetTitleOffset(x_title_offset); DT_histo2->GetYaxis()->SetTitleOffset(y_title_offset);   DT_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    DT_histo2->GetXaxis()->SetTitle(xaxislabel.c_str());   DT_histo2->GetYaxis()->SetTitle(yaxislabel.c_str());     DT_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "DT");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c4, pdfname4.c_str());

    y_range_2 = 10.00;
    std::stringstream pdfname5ss; pdfname5ss<<"STAMuonHits_2D_2023_"<<"RedPercent"<<"_CSC"; std::string pdfname5 = pdfname5ss.str(); 
    TCanvas * c5 = new TCanvas(pdfname5.c_str(), "", 800, 400); c5->cd(); c5->SetTicks(1,1); c5->SetLeftMargin(0.05); c5->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    CSC_histo2->SetMarkerSize(1.75); CSC_histo2->Draw("textcolz"); CSC_histo2->GetYaxis()->SetRangeUser(y_range_1,y_range_2); CSC_histo2->GetZaxis()->SetRangeUser(0,100);
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "CSC");
    CSC_histo2->GetXaxis()->SetTitleOffset(x_title_offset); CSC_histo2->GetYaxis()->SetTitleOffset(y_title_offset);   CSC_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    CSC_histo2->GetXaxis()->SetTitle(xaxislabel.c_str());   CSC_histo2->GetYaxis()->SetTitle(yaxislabel.c_str());     CSC_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c5, pdfname5.c_str());
    */

    double pogbinning[] = { -2.4,  -2.1, -1.6,  -1.2,  -0.9,  -0.6,  -0.3,  -0.2,   0.2,   0.3,   0.6,   0.9,   1.2,   1.6,   2.1,   2.4};
    double tightideff[] = {    0.933, 0.924, 0.969, 0.960, 0.963, 0.977, 0.876, 0.962, 0.868, 0.976, 0.961, 0.958, 0.960, 0.944, 0.944  };

    TH1F * TightIDEff = new TH1F("TightIDEff", "2012 Tight ID Efficiency", 15, pogbinning);
    for(int i=0; i<15; ++i) { TightIDEff->SetBinContent(i+1, tightideff[i]); }

    /*
    {'-0.2_0.2': {'data': {'efficiency': 0.9620470334987958,
	'-0.3_-0.2': {'data': {'efficiency': 0.8785896650095353,
	'-0.6_-0.3': {'data': {'efficiency': 0.9769343813246657,
	'-0.9_-0.6': {'data': {'efficiency': 0.9630945761265595,
	'-1.2_-0.9': {'data': {'efficiency': 0.9603076465268146,
	'-1.6_-1.2': {'data': {'efficiency': 0.969156254195653,
	'-2.1_-1.6': {'data': {'efficiency': 0.9236438999539651,
	'-2.4_-2.1': {'data': {'efficiency': 0.9329682595251612,
	'0.2_0.3': {'data': {'efficiency': 0.8681874065897149,
	'0.3_0.6': {'data': {'efficiency': 0.9757265466689626,
	'0.6_0.9': {'data': {'efficiency': 0.9614116022156562,
	'0.9_1.2': {'data': {'efficiency': 0.9575638243542479,
	'1.2_1.6': {'data': {'efficiency': 0.9605016006018594,
	'1.6_2.1': {'data': {'efficiency': 0.9442083234514784,
	'2.1_2.4': {'data': {'efficiency': 0.9437111925674662,
    */


    // reduced eta binning
    // double reduced[] = {-2.5, -1.7, -1.6, -1.1, -0.9, -0.8, -0.3, -0.2, 0.2, 0.3, 0.8, 0.9, 1.1, 1.6, 1.7, 2.5};
    double reduced[] = { -2.4,  -2.1, -1.6,  -1.2,  -0.9,  -0.6,  -0.3,  -0.2,   0.2,   0.3,   0.6,   0.9,   1.2,   1.6,   2.1,   2.4};

    // TH2F * Efficiency_2D_Plot_100 = new TH2F("Efficiency_2D_Plot_100", "All Hits/Segments 100% Efficient", 15, reduced, 16, -0.5, 15.5);
    // TH2F * Efficiency_2D_Plot_095 = new TH2F("Efficiency_2D_Plot_095", "All Hits/Segments  95% Efficient", 15, reduced, 16, -0.5, 15.5);
    // TH2F * Efficiency_2D_Plot_090 = new TH2F("Efficiency_2D_Plot_090", "All Hits/Segments  90% Efficient", 15, reduced, 16, -0.5, 15.5);
    // TH2F * Efficiency_2D_Plot_080 = new TH2F("Efficiency_2D_Plot_080", "All Hits/Segments  80% Efficient", 15, reduced, 16, -0.5, 15.5);
    // TH2F * Efficiency_2D_Plot_050 = new TH2F("Efficiency_2D_Plot_050", "All Hits/Segments  50% Efficient", 15, reduced, 16, -0.5, 15.5);
    // std::vector< TH2F * > Efficiency_2D_Plot; 
    // Efficiency_2D_Plot.push_back(Efficiency_2D_Plot_100); Efficiency_2D_Plot.push_back(Efficiency_2D_Plot_095); Efficiency_2D_Plot.push_back(Efficiency_2D_Plot_090);
    // Efficiency_2D_Plot.push_back(Efficiency_2D_Plot_080); Efficiency_2D_Plot.push_back(Efficiency_2D_Plot_050);

    TH1F * Efficiency_1D_Plot_100 = new TH1F("Efficiency_1D_Plot_100", "All Hits/Segments 100% Efficient", 15, reduced);
    TH1F * Efficiency_1D_Plot_095 = new TH1F("Efficiency_1D_Plot_095", "All Hits/Segments  95% Efficient", 15, reduced);
    TH1F * Efficiency_1D_Plot_090 = new TH1F("Efficiency_1D_Plot_090", "All Hits/Segments  90% Efficient", 15, reduced);
    TH1F * Efficiency_1D_Plot_080 = new TH1F("Efficiency_1D_Plot_080", "All Hits/Segments  80% Efficient", 15, reduced);
    TH1F * Efficiency_1D_Plot_050 = new TH1F("Efficiency_1D_Plot_050", "All Hits/Segments  50% Efficient", 15, reduced);
    TH1F * Efficiency_1D_Plot_SS1 = new TH1F("Efficiency_1D_Plot_SS1", "All Hits/Segments Special Scenario 1", 15, reduced);
    TH1F * Efficiency_1D_Plot_SS2 = new TH1F("Efficiency_1D_Plot_SS2", "All Hits/Segments Special Scenario 2", 15, reduced);
    std::vector< TH1F * > Efficiency_1D_Plot; 
    Efficiency_1D_Plot.push_back(Efficiency_1D_Plot_100); Efficiency_1D_Plot.push_back(Efficiency_1D_Plot_095); Efficiency_1D_Plot.push_back(Efficiency_1D_Plot_090);
    Efficiency_1D_Plot.push_back(Efficiency_1D_Plot_080); Efficiency_1D_Plot.push_back(Efficiency_1D_Plot_050); 
    Efficiency_1D_Plot.push_back(Efficiency_1D_Plot_SS1); Efficiency_1D_Plot.push_back(Efficiency_1D_Plot_SS2);


    // segment / rechit efficiency
    // double hiteff[] = {1.00, 0.95, 0.90, 0.80, 0.50};
    double hiteff[] = {0.90, 0.80, 0.70, 0.60, 0.50, 0.90, 0.80};


    double two_out_of_X_eff[7][15]; 

    for(int k=0; k<5; ++k) {
      for(int j=2; j<17; ++j) { 

	// Eff = 1.00 - Probability only 1 segment is reconstructed                - Probability zero segments are reconstructed
        //     = 1.00 - Combination of 1 out of j x hiteff^1 x (1 - hiteff)^(j-1)  - Combination of 0 out of j x hiteff^0 x (1 - hiteff)^(j-0)
	// 
	// Full Formula:
	// -------------
	// double eff = 1 
	//   - factorial(j) * 1.0 / (factorial(1) * factorial(j-1)) * pow(hiteff[k],1) * pow(1-hiteff[k],j-1)
	//   - factorial(j) * 1.0 / (factorial(0) * factorial(j-0)) * pow(hiteff[k],0) * pow(1-hiteff[k],j-0);
	//
	// Reduced Formula:
	// ----------------
	double eff = 1 - factorial(j) * 1.0 / factorial(j-1) * hiteff[k] * pow(1-hiteff[k],j-1) - 1.0 * pow(1-hiteff[k],j);

	std::cout<<"Loop ["<<k<<"]["<<j<<"] :: Amount of hits = "<<j<<" efficiency of single hit = "<<hiteff[k]<<" results in total efficiency = "<<eff<<std::endl;
	two_out_of_X_eff[k][j-2] = eff;
      }
    }

    // Special Scenario 1 :: all chambers at 90% efficiency + 1 chamber at 0% 
    // Special Scenario 2 :: all chambers at 80% efficiency + 1 chamber at 0%
    // Strategy :: do j-1; you automatically loose one station that does not enter into the efficiency calculations;
    //             therefore start with j = 3, because j == 2 is automatically inefficient
    two_out_of_X_eff[5][0] = 0; two_out_of_X_eff[6][0] = 0;
    //             start loop at j = 3 and do the calculations for (j-1)
    for(int k=5; k<7; ++k) {
      for(int j=3; j<17; ++j) {
	double eff = 1 - factorial(j-1) * 1.0 / factorial(j-2) * hiteff[k] * pow(1-hiteff[k],j-2) - 1.0 * pow(1-hiteff[k],j-1);
	std::cout<<"Loop ["<<k<<"]["<<j<<"] :: Inactive hit = 1, Active hits = "<<j-1<<" efficiency of single hit = "<<hiteff[k]<<" results in total efficiency = "<<eff<<std::endl;
	two_out_of_X_eff[k][j-2] = eff;
      }
    }

    for(int k=0; k<7; ++k) { 
      for(int i=0; i<15; ++i) {
	double final_eff = 0.0;
	for(int  j=2; j<17; ++j) {
	  double weight = All_histo2->GetBinContent(i+1,j+1);
	  final_eff += weight * 1.0/100.0 * two_out_of_X_eff[k][j-2];
	  // std::cout<<"Loop ["<<k<<"]["<<i<<"]["<<j<<"] :: Amount of hits = "<<j<<" with weight = "<<weight<<" % and efficiency of single hit = "<<hiteff[k];
	  // std::cout<<" efficiency of 2 out of "<<j<<" = "<<two_out_of_X_eff[k][j-2]<<" results in total efficiency = "<<final_eff<<std::endl;
	}
	Efficiency_1D_Plot[k]->SetBinContent(i+1, final_eff*TightIDEff->GetBinContent(i+1));
      }
    }

    TightIDEff->SetLineColor(kBlack);               TightIDEff->SetMarkerColor(kBlack);
    /* Original Colors for showing 90% 80% 70% 60% 50% */ 
    /*
    Efficiency_1D_Plot[0]->SetLineColor(kGreen+2);  Efficiency_1D_Plot[0]->SetMarkerColor(kGreen+2);
    Efficiency_1D_Plot[1]->SetLineColor(kOrange);   Efficiency_1D_Plot[1]->SetMarkerColor(kOrange);
    Efficiency_1D_Plot[2]->SetLineColor(kBlue);     Efficiency_1D_Plot[2]->SetMarkerColor(kBlue);
    Efficiency_1D_Plot[3]->SetLineColor(kViolet);   Efficiency_1D_Plot[3]->SetMarkerColor(kViolet);
    Efficiency_1D_Plot[4]->SetLineColor(kRed);      Efficiency_1D_Plot[4]->SetMarkerColor(kRed);
    */
    Efficiency_1D_Plot[0]->SetLineColor(kGreen+2);  Efficiency_1D_Plot[0]->SetMarkerColor(kGreen+2);    // 90%
    Efficiency_1D_Plot[1]->SetLineColor(kBlue);     Efficiency_1D_Plot[1]->SetMarkerColor(kBlue);       // 80% 
    // Efficiency_1D_Plot[2]->SetLineColor(kBlue);     Efficiency_1D_Plot[2]->SetMarkerColor(kBlue);    // 70%
    Efficiency_1D_Plot[3]->SetLineColor(kRed);      Efficiency_1D_Plot[3]->SetMarkerColor(kRed);        // 60%
    // Efficiency_1D_Plot[4]->SetLineColor(kRed);      Efficiency_1D_Plot[4]->SetMarkerColor(kRed);     // 50%

    TightIDEff->SetLineWidth(2); TightIDEff->SetMarkerStyle(34); TightIDEff->SetMarkerSize(1.50);
    for(int k=0; k<5; ++k) { Efficiency_1D_Plot[k]->SetLineWidth(2); Efficiency_1D_Plot[k]->SetMarkerStyle(20); Efficiency_1D_Plot[k]->SetMarkerSize(1.00); }

    // New: different markerstyles for 80% and 60%
    // circle:  20 & 24
    // square:  21 & 25
    // diamond: 33 & 27 --> rediculously small --> markersize 1.5
    Efficiency_1D_Plot[1]->SetMarkerStyle(21);
    Efficiency_1D_Plot[3]->SetMarkerStyle(33); Efficiency_1D_Plot[3]->SetMarkerSize(1.50);

    /* Original Colors for showing 90% 80% 70% 60% 50% */ 
    /*
    Efficiency_1D_Plot[5]->SetLineColor(kGreen+2);  Efficiency_1D_Plot[5]->SetMarkerColor(kGreen+2); Efficiency_1D_Plot[5]->SetLineWidth(2); Efficiency_1D_Plot[5]->SetLineStyle(2); Efficiency_1D_Plot[5]->SetMarkerStyle(24);
    Efficiency_1D_Plot[6]->SetLineColor(kOrange);   Efficiency_1D_Plot[6]->SetMarkerColor(kOrange);  Efficiency_1D_Plot[6]->SetLineWidth(2); Efficiency_1D_Plot[6]->SetLineStyle(2); Efficiency_1D_Plot[6]->SetMarkerStyle(24);
    */
    Efficiency_1D_Plot[5]->SetLineColor(kGreen+2);  Efficiency_1D_Plot[5]->SetMarkerColor(kGreen+2); Efficiency_1D_Plot[5]->SetLineWidth(2); Efficiency_1D_Plot[5]->SetLineStyle(2); Efficiency_1D_Plot[5]->SetMarkerStyle(24);
    Efficiency_1D_Plot[6]->SetLineColor(kBlue);     Efficiency_1D_Plot[6]->SetMarkerColor(kBlue);    Efficiency_1D_Plot[6]->SetLineWidth(2); Efficiency_1D_Plot[6]->SetLineStyle(2); Efficiency_1D_Plot[6]->SetMarkerStyle(24);
    // New: different markerstyles for 80% and 60%
    Efficiency_1D_Plot[6]->SetMarkerStyle(25);

    y_range_1 =  0.60; // 0.40 Originally
    y_range_2 =  1.00;
    std::string yaxislabel_ineff = "RECO + ID Efficiency";
    std::string xaxislabel_ineff = "#eta";
    double x_title_offset_ineff = 0.75;
    double y_title_offset_ineff = 1.00;
    // std::string info5 = "Stand Alone Muon Reco Efficiency";
    // std::string info6 = "one RPC hit / DT,CSC segment missing";
    y_title_offset = 1.25;
    // Ineff_All_histo->SetLineWidth(2.0); Ineff_All_histo->SetMarkerStyle(20); Ineff_All_histo->SetMarkerSize(1.25);
    std::stringstream pdfname6ss; pdfname6ss<<"STAMuonHits_2D_2023_"<<"GradualEfficiency"<<"_All"; std::string pdfname6 = pdfname6ss.str();
    TCanvas * c6 = new TCanvas(pdfname6.c_str(), "", 800, 400); c6->cd(); c6->SetTicks(1,1); c6->SetLeftMargin(0.075); c6->SetRightMargin(0.025); // gStyle->SetPaintTextFormat("3.0f");
    TightIDEff->Draw("HP");
    Efficiency_1D_Plot[0]->Draw("HPSame"); Efficiency_1D_Plot[1]->Draw("HPSame"); /*Efficiency_1D_Plot[2]->Draw("HPSame");*/ Efficiency_1D_Plot[3]->Draw("HPSame"); /*Efficiency_1D_Plot[4]->Draw("HPSame"); */
    latex_cmslab.DrawLatex(0.075, 0.925,"CMS Phase 2 Simulation"); /*#sqrt{s} = 13 TeV*/ /*latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str());*/ gPad->RedrawAxis(); // latex_cmslab.DrawLatex(0.125, 0.85, "All");
    TightIDEff->GetXaxis()->SetTitleOffset(x_title_offset_ineff); TightIDEff->GetYaxis()->SetTitleOffset(y_title_offset_ineff);      
    TightIDEff->GetXaxis()->SetTitle(xaxislabel_ineff.c_str());   TightIDEff->GetYaxis()->SetTitle(yaxislabel_ineff.c_str());        
    TightIDEff->GetXaxis()->SetTitleSize(0.06);
    TightIDEff->GetYaxis()->SetRangeUser(y_range_1,y_range_2);
    c6->Update();
    TightIDEff->GetYaxis()->SetNdivisions(504, kFALSE);
    c6->Update();
    gPad->RedrawAxis();
    c6->Update();
    // Originally:     double leg2_x1 = 0.30, leg2_x2 = 0.70, leg2_y1 = 0.15, leg2_y2 = 0.50;
    double leg2_x1 = 0.30, leg2_x2 = 0.70, leg2_y1 = 0.15, leg2_y2 = 0.40;
    TLegend * legend = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2, NULL,"brNDC"); legend->SetLineColor(1);    legend->SetLineStyle(1);  legend->SetLineWidth(1);  legend->SetFillColor(4000); legend->SetBorderSize(1);
    legend->AddEntry(TightIDEff,            "Tight ID (2012) Efficiency",   "fp");
    legend->AddEntry(Efficiency_1D_Plot[0], "Segments/hits  90% Efficient", "fp");
    legend->AddEntry(Efficiency_1D_Plot[1], "Segments/hits  80% Efficient", "fp");
    // legend->AddEntry(Efficiency_1D_Plot[2], "Segments/hits  70% Efficient", "fp");
    legend->AddEntry(Efficiency_1D_Plot[3], "Segments/hits  60% Efficient", "fp");
    // legend->AddEntry(Efficiency_1D_Plot[4], "Segments/hits  50% Efficient", "fp");
    legend->Draw();
    c6->Update();
    savePlot(c6, pdfname6.c_str());

    std::stringstream pdfname7ss; pdfname7ss<<"STAMuonHits_2D_2023_"<<"GradualEfficiency"<<"_All_v2"; std::string pdfname7 = pdfname7ss.str();
    TCanvas * c7 = new TCanvas(pdfname7.c_str(), "", 800, 600); c7->cd(); c7->SetTicks(1,1); // c7->SetLeftMargin(0.05); c7->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    TightIDEff->Draw("HP");
    Efficiency_1D_Plot[0]->Draw("HPSame"); Efficiency_1D_Plot[1]->Draw("HPSame"); /*Efficiency_1D_Plot[2]->Draw("HPSame");*/ Efficiency_1D_Plot[3]->Draw("HPSame"); /*Efficiency_1D_Plot[4]->Draw("HPSame");*/ 
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Phase 2 Simulation"); /*#sqrt{s} = 13 TeV*/ /*latex_cmslab.DrawLatex(0.575, 0.925, label2.c_str());*/ gPad->RedrawAxis(); // latex_cmslab.DrawLatex(0.125, 0.85, "All");
    TightIDEff->GetXaxis()->SetTitleOffset(x_title_offset_ineff); TightIDEff->GetYaxis()->SetTitleOffset(y_title_offset_ineff);
    TightIDEff->GetXaxis()->SetTitle(xaxislabel_ineff.c_str());   TightIDEff->GetYaxis()->SetTitle(yaxislabel_ineff.c_str());
    TightIDEff->GetXaxis()->SetTitleSize(0.06);
    TightIDEff->GetYaxis()->SetRangeUser(y_range_1,y_range_2);
    TightIDEff->GetYaxis()->SetNdivisions(504, kFALSE);
    c7->Update();
    legend->Draw();
    c7->Update();
    savePlot(c7, pdfname7.c_str());

    std::stringstream pdfname9ss; pdfname9ss<<"STAMuonHits_2D_2023_"<<"GradualEfficiency"<<"_All"; std::string pdfname9 = pdfname9ss.str();
    TCanvas * c9 = new TCanvas(pdfname9.c_str(), "", 800, 400); c9->cd(); c9->SetTicks(1,1); c9->SetLeftMargin(0.075); c9->SetRightMargin(0.025); // gStyle->SetPaintTextFormat("3.0f");
    TightIDEff->Draw("HP");
    Efficiency_1D_Plot[0]->Draw("HPSame"); Efficiency_1D_Plot[1]->Draw("HPSame"); /*Efficiency_1D_Plot[2]->Draw("HPSame");*/ Efficiency_1D_Plot[3]->Draw("HPSame"); /*Efficiency_1D_Plot[4]->Draw("HPSame"); */
    latex_cmslab.DrawLatex(0.075, 0.925,"CMS Phase 2 Simulation"); /*#sqrt{s} = 13 TeV*/ /*latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str());*/ gPad->RedrawAxis(); // latex_cmslab.DrawLatex(0.125, 0.85, "All");
    TightIDEff->GetXaxis()->SetTitleOffset(x_title_offset_ineff); TightIDEff->GetYaxis()->SetTitleOffset(y_title_offset_ineff);      
    TightIDEff->GetXaxis()->SetTitle(xaxislabel_ineff.c_str());   TightIDEff->GetYaxis()->SetTitle(yaxislabel_ineff.c_str());        
    TightIDEff->GetXaxis()->SetTitleSize(0.06);
    TightIDEff->GetYaxis()->SetRangeUser(y_range_1,y_range_2);
    TightIDEff->GetYaxis()->SetNdivisions(504, kFALSE);
    c9->Update();
    legend->Draw();
    c9->Update();
    savePlot(c9, pdfname9.c_str());

    std::stringstream pdfname8ss; pdfname8ss<<"STAMuonHits_2D_2023_"<<"GradualEfficiency"<<"_All_SpecialScenarios"; std::string pdfname8 = pdfname8ss.str();
    TCanvas * c8 = new TCanvas(pdfname8.c_str(), "", 800, 400); c8->cd(); c8->SetTicks(1,1); c8->SetLeftMargin(0.075); c8->SetRightMargin(0.025); // gStyle->SetPaintTextFormat("3.0f");
    TightIDEff->Draw("HP");
    Efficiency_1D_Plot[0]->Draw("HPSame"); Efficiency_1D_Plot[1]->Draw("HPSame"); /*Efficiency_1D_Plot[2]->Draw("HPSame");*/ Efficiency_1D_Plot[3]->Draw("HPSame"); /*Efficiency_1D_Plot[4]->Draw("HPSame");*/ 
    Efficiency_1D_Plot[5]->Draw("HPSame"); Efficiency_1D_Plot[6]->Draw("HPSame");
    latex_cmslab.DrawLatex(0.075, 0.925,"CMS Phase 2 Simulation"); /*#sqrt{s} = 13 TeV*/ /*latex_cmslab.DrawLatex(0.70, 0.925, label2.c_str()); gPad->RedrawAxis();*/ // latex_cmslab.DrawLatex(0.125, 0.85, "All");
    TightIDEff->GetXaxis()->SetTitleOffset(x_title_offset_ineff); TightIDEff->GetYaxis()->SetTitleOffset(y_title_offset_ineff);
    TightIDEff->GetXaxis()->SetTitle(xaxislabel_ineff.c_str());   TightIDEff->GetYaxis()->SetTitle(yaxislabel_ineff.c_str());
    TightIDEff->GetXaxis()->SetTitleSize(0.06);
    TightIDEff->GetYaxis()->SetRangeUser(y_range_1,y_range_2);
    TightIDEff->GetYaxis()->SetNdivisions(504, kFALSE);
    c8->Update();
    legend->AddEntry(Efficiency_1D_Plot[5], "1 Station off, Others 90% Efficient", "fp");
    legend->AddEntry(Efficiency_1D_Plot[6], "1 Station off, Others 80% Efficient", "fp");
    legend->Draw();
    c8->Update();
    savePlot(c8, pdfname8.c_str());


}

int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void savePlot(TCanvas * c, TString name) {
  c->SaveAs(name+".png");
  c->SaveAs(name+".eps");
  gSystem->Exec("epstopdf "+name+".eps");
  std::cout<<"Plot "<<name<<" saved"<<std::endl;
}
