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

void Percentage_2D_StaMuonHits_2023() {

    // Title
    // -----
    std::string title = "62Xvs71X";


    // Main Comparison
    // ---------------
    std::string label1 = "CMSSW 62X SLHC (2023 Geometry)"; 
    std::string label2 = "CMSSW 71X (2023 Geometry)"; 

    // Additional Information
    // ----------------------
    std::string info1 = "Hits used in Stand Alone Muon Fit";
    std::string info2 = "these are Segments (CSC,DT) and Rechits (RPC)";
    std::string info3 = "Sample used: RelValZmumuJets_Pt_20_300";

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
    std::string Histo_RPCb = "SegmentsAndHits/Rechits_All_Eta_2D";
    std::string Histo_RPCf = "SegmentsAndHits/Rechits_RPC_Eta_2D";
    std::string Histo_CSC  = "SegmentsAndHits/Rechits_CSC_Eta_2D";
    std::string Histo_DT   = "SegmentsAndHits/Rechits_DT_Eta_2D";

    std::string Histo_All  = "Muon_All";

    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1,0);
    gStyle->SetOptStat(000000);
    gStyle->SetOptTitle(0);
    gStyle->SetGridStyle(1);    // grid
    gROOT->UseCurrentStyle();
    
    TFile * infile1             = new TFile((RootFileName1).c_str(),"READ");     std::cout<<"File 1  = "<<infile1<<std::endl; 
    TFile * infile2             = new TFile((RootFileName2).c_str(),"READ");     std::cout<<"File 2  = "<<infile2<<std::endl; 

    TH2F * RPCb_histo1          = (TH2F*) infile1->Get(Histo_RPCb.c_str());  std::cout<<"Histo  = "<<RPCb_histo1<<std::endl; 
    TH2F * RPCb_histo2          = (TH2F*) infile2->Get(Histo_RPCb.c_str());  std::cout<<"Histo  = "<<RPCb_histo2<<std::endl; 
    TH2F * RPCf_histo1          = (TH2F*) infile1->Get(Histo_RPCf.c_str());  std::cout<<"Histo  = "<<RPCf_histo1<<std::endl; 
    TH2F * RPCf_histo2          = (TH2F*) infile2->Get(Histo_RPCf.c_str());  std::cout<<"Histo  = "<<RPCf_histo2<<std::endl; 
    TH2F * CSC_histo1           = (TH2F*) infile1->Get(Histo_CSC.c_str());   std::cout<<"Histo  = "<<CSC_histo1<<std::endl; 
    TH2F * CSC_histo2           = (TH2F*) infile2->Get(Histo_CSC.c_str());   std::cout<<"Histo  = "<<CSC_histo2<<std::endl; 
    TH2F * DT_histo1            = (TH2F*) infile1->Get(Histo_DT.c_str());    std::cout<<"Histo  = "<<DT_histo1<<std::endl; 
    TH2F * DT_histo2            = (TH2F*) infile2->Get(Histo_DT.c_str());    std::cout<<"Histo  = "<<DT_histo2<<std::endl; 

    TH1F * All_histo1            = (TH1F*) infile1->Get(Histo_All.c_str());    std::cout<<"Histo  = "<<All_histo1<<std::endl; 
    TH1F * All_histo2            = (TH1F*) infile2->Get(Histo_All.c_str());    std::cout<<"Histo  = "<<All_histo2<<std::endl; 

    RPCb_histo1->SetLineColor(kGray+2); RPCb_histo1->SetLineWidth(1);         RPCb_histo2->SetLineColor(kRed); RPCb_histo2->SetLineWidth(1); 
    RPCb_histo1->SetFillColor(kBlue-7); /*RPCb_histo1->SetFillStyle(4050);*/  RPCb_histo2->SetFillColor(kRed); RPCb_histo2->SetFillStyle(3002); 
    RPCf_histo1->SetLineColor(kGray+2); RPCf_histo1->SetLineWidth(1);         RPCf_histo2->SetLineColor(kRed); RPCf_histo2->SetLineWidth(1); 
    RPCf_histo1->SetFillColor(kBlue-7); /*RPCf_histo1->SetFillStyle(4050);*/  RPCf_histo2->SetFillColor(kRed); RPCf_histo2->SetFillStyle(3002); 
    CSC_histo1->SetLineColor(kGray+2);  CSC_histo1->SetLineWidth(1);          CSC_histo2->SetLineColor(kRed);  CSC_histo2->SetLineWidth(1); 
    CSC_histo1->SetFillColor(kBlue-7);  /*CSC_histo1->SetFillStyle(4050);*/   CSC_histo2->SetFillColor(kRed);  CSC_histo2->SetFillStyle(3002); 
    DT_histo1->SetLineColor(kGray+2);   DT_histo1->SetLineWidth(1);           DT_histo2->SetLineColor(kRed);   DT_histo2->SetLineWidth(1); 
    DT_histo1->SetFillColor(kBlue-7);   /*DT_histo1->SetFillStyle(4050);*/    DT_histo2->SetFillColor(kRed);   DT_histo2->SetFillStyle(3002); 
    // std::string Xaxis = "log_{10} TOF (ns)";
    std::string Xaxis = "#eta";
    std::string Yaxis = "\# Hits";

    RPCb_histo1->SetTitle("All");  RPCb_histo1->GetXaxis()->SetTitle(Xaxis.c_str()); RPCb_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    RPCf_histo1->SetTitle("RPC");  RPCf_histo1->GetXaxis()->SetTitle(Xaxis.c_str()); RPCf_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    CSC_histo1->SetTitle("CSC");   CSC_histo1->GetXaxis()->SetTitle(Xaxis.c_str());  CSC_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    DT_histo1->SetTitle("DT");     DT_histo1->GetXaxis()->SetTitle(Xaxis.c_str());   DT_histo1->GetYaxis()->SetTitle(Yaxis.c_str());

    double denominator, numerator;
    for(int i=0; i<RPCb_histo1->GetXaxis()->GetNbins(); ++i) {
      denominator = All_histo1->GetBinContent(i+1);
      for(int j=0; j<RPCb_histo1->GetYaxis()->GetNbins(); ++j) {
	numerator = RPCb_histo1->GetBinContent(i+1,j+1);
	if (numerator > 0.0 && denominator > 0.0) RPCb_histo1->SetBinContent(i+1,j+1, numerator/denominator*100);
	numerator = RPCf_histo1->GetBinContent(i+1,j+1);
	if (numerator > 0.0 && denominator > 0.0) RPCf_histo1->SetBinContent(i+1,j+1, numerator/denominator*100);
	numerator = CSC_histo1->GetBinContent(i+1,j+1);
	if (numerator > 0.0 && denominator > 0.0) CSC_histo1->SetBinContent(i+1,j+1, numerator/denominator*100);
	numerator = DT_histo1->GetBinContent(i+1,j+1);
	if (numerator > 0.0 && denominator > 0.0) DT_histo1->SetBinContent(i+1,j+1, numerator/denominator*100);
      }
    }


    int RPCb_1_entries = RPCb_histo1->GetEntries();       int RPCb_2_entries = RPCb_histo2->GetEntries();
    int RPCf_1_entries = RPCf_histo1->GetEntries();       int RPCf_2_entries = RPCf_histo2->GetEntries();
    int CSC_1_entries  = CSC_histo1->GetEntries();        int CSC_2_entries  = CSC_histo2->GetEntries();
    int DT_1_entries   = DT_histo1->GetEntries();         int DT_2_entries   = DT_histo2->GetEntries();

    std::stringstream ss_RPC_b_1_label; ss_RPC_b_1_label << label1;      std::stringstream ss_RPC_b_2_label; ss_RPC_b_2_label << label2;
    std::stringstream ss_RPC_f_1_label; ss_RPC_f_1_label << label1;      std::stringstream ss_RPC_f_2_label; ss_RPC_f_2_label << label2;
    std::stringstream ss_CSC_1_label;   ss_CSC_1_label   << label1;      std::stringstream ss_CSC_2_label;   ss_CSC_2_label   << label2;
    std::stringstream ss_DT_1_label;    ss_DT_1_label    << label1;      std::stringstream ss_DT_2_label;    ss_DT_2_label    << label2;

    std::string RPC_b_1_label = ss_RPC_b_1_label.str();    std::string RPC_b_2_label = ss_RPC_b_2_label.str();
    std::string RPC_f_1_label = ss_RPC_f_1_label.str();    std::string RPC_f_2_label = ss_RPC_f_2_label.str();
    std::string CSC_1_label   = ss_CSC_1_label.str();      std::string CSC_2_label = ss_CSC_2_label.str();
    std::string DT_1_label    = ss_DT_1_label.str();       std::string DT_2_label = ss_DT_2_label.str();

    /*
    RPCb_histo1->Scale(1.0/RPCb_1_entries); RPCb_histo2->Scale(1.0/RPCb_2_entries);
    RPCf_histo1->Scale(1.0/RPCf_1_entries); RPCf_histo2->Scale(1.0/RPCf_2_entries);
    CSC_histo1->Scale(1.0/CSC_1_entries);   CSC_histo2->Scale(1.0/CSC_2_entries);
    DT_histo1->Scale(1.0/DT_1_entries);     DT_histo2->Scale(1.0/DT_2_entries);
    */

    double leg_x1 = 0.55, leg_x2 = 0.85, leg_y1 = 0.700, leg_y2 = 0.850; 

    TLegend * l1 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l1->SetLineColor(1);    l1->SetLineStyle(1);  l1->SetLineWidth(1);  l1->SetFillColor(4000); l1->SetBorderSize(1);
    TLegend * l2 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l2->SetLineColor(1);    l2->SetLineStyle(1);  l2->SetLineWidth(1);  l2->SetFillColor(4000); l2->SetBorderSize(1);
    TLegend * l3 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l3->SetLineColor(1);    l3->SetLineStyle(1);  l3->SetLineWidth(1);  l3->SetFillColor(4000); l3->SetBorderSize(1);
    TLegend * l4 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l4->SetLineColor(1);    l4->SetLineStyle(1);  l4->SetLineWidth(1);  l4->SetFillColor(4000); l4->SetBorderSize(1);

    l1->AddEntry(RPCb_histo1, RPC_b_1_label.c_str(), "F");     l1->AddEntry(RPCb_histo2, RPC_b_2_label.c_str(), "F");
    l2->AddEntry(RPCf_histo1, RPC_f_1_label.c_str(), "F");     l2->AddEntry(RPCf_histo2, RPC_f_2_label.c_str(), "F");
    l3->AddEntry(CSC_histo1,  CSC_1_label.c_str(),   "F");     l3->AddEntry(CSC_histo2,  CSC_2_label.c_str(),   "F");
    l4->AddEntry(DT_histo1,   DT_1_label.c_str(),    "F");     l4->AddEntry(DT_histo2,   DT_2_label.c_str(),    "F");

    TLatex latex_cmslab; latex_cmslab.SetNDC(); latex_cmslab.SetTextSize(0.035); latex_cmslab.SetTextAlign(11);
    TLatex latex_info;   latex_info.SetNDC();   latex_info.SetTextSize(0.025);    latex_info.SetTextAlign(31);

    double infoX_x_right = 0.6275; 
    double info1_y = 0.80, info2_y = 0.75, info3_y = 0.70;

    /*
    TLegend * l5 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l5->SetLineColor(1);    l5->SetLineStyle(1);  l5->SetLineWidth(1);  l5->SetFillColor(4000); l5->SetBorderSize(1);
    TLegend * l6 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l6->SetLineColor(1);    l6->SetLineStyle(1);  l6->SetLineWidth(1);  l6->SetFillColor(4000); l6->SetBorderSize(1);
    TLegend * l7 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l7->SetLineColor(1);    l7->SetLineStyle(1);  l7->SetLineWidth(1);  l7->SetFillColor(4000); l7->SetBorderSize(1);
    TLegend * l8 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l8->SetLineColor(1);    l8->SetLineStyle(1);  l8->SetLineWidth(1);  l8->SetFillColor(4000); l8->SetBorderSize(1);

    l5->AddEntry(RPCb_histo1, RPC_b_1_label.c_str(), "F");     l5->AddEntry(RPCb_histo2, RPC_b_2_label.c_str(), "F");
    l6->AddEntry(RPCf_histo1, RPC_f_1_label.c_str(), "F");     l6->AddEntry(RPCf_histo2, RPC_f_2_label.c_str(), "F");
    l7->AddEntry(CSC_histo1,  CSC_1_label.c_str(),   "F");     l7->AddEntry(CSC_histo2,  CSC_2_label.c_str(),   "F");
    l8->AddEntry(DT_histo1,   DT_1_label.c_str(),    "F");     l8->AddEntry(DT_histo2,   DT_2_label.c_str(),    "F");
    */

    // Save the Other Canvasses
    y_range_2 = 15.00;
    std::stringstream pdfname2ss; pdfname2ss<<"STAMuonHits_2D_2023_"<<"Percentage"<<"_All"; std::string pdfname2 = pdfname2ss.str(); 
    TCanvas * c2 = new TCanvas(pdfname2.c_str(), "", 800, 400); c2->cd(); c2->SetTicks(1,1); c2->SetLeftMargin(0.05); c2->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    RPCb_histo1->Draw("textcolz"); /*l5->Draw();*/ RPCb_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); RPCb_histo1->GetZaxis()->SetRangeUser(0,100); 
    RPCb_histo1->GetXaxis()->SetTitleOffset(x_title_offset); RPCb_histo1->GetYaxis()->SetTitleOffset(y_title_offset);   RPCb_histo1->GetZaxis()->SetTitleOffset(z_title_offset);
    RPCb_histo1->GetXaxis()->SetTitle(xaxislabel.c_str());   RPCb_histo1->GetYaxis()->SetTitle(yaxislabel.c_str());     RPCb_histo1->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label1.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "All");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c2, pdfname2.c_str());

    y_range_2 = 10.00;
    std::stringstream pdfname3ss; pdfname3ss<<"STAMuonHits_2D_2023_"<<"Percentage"<<"_RPC"; std::string pdfname3 = pdfname3ss.str(); 
    TCanvas * c3 = new TCanvas(pdfname3.c_str(), "", 800, 400); c3->cd(); c3->SetTicks(1,1); c3->SetLeftMargin(0.05); c3->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    RPCf_histo1->Draw("textcolz"); /*l6->Draw();*/ RPCf_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); RPCf_histo1->GetZaxis()->SetRangeUser(0,100);
    RPCf_histo1->GetXaxis()->SetTitleOffset(x_title_offset); RPCf_histo1->GetYaxis()->SetTitleOffset(y_title_offset);   RPCf_histo1->GetZaxis()->SetTitleOffset(z_title_offset);
    RPCf_histo1->GetXaxis()->SetTitle(xaxislabel.c_str());   RPCf_histo1->GetYaxis()->SetTitle(yaxislabel.c_str());     RPCf_histo1->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label1.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "RPC");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c3, pdfname3.c_str());

    y_range_2 = 10.00;
    std::stringstream pdfname4ss; pdfname4ss<<"STAMuonHits_2D_2023_"<<"Percentage"<<"_DT"; std::string pdfname4 = pdfname4ss.str(); 
    TCanvas * c4 = new TCanvas(pdfname4.c_str(), "", 800, 400); c4->cd(); c4->SetTicks(1,1); c4->SetLeftMargin(0.05); c4->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    DT_histo1->Draw("textcolz"); /*l8->Draw();*/ DT_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); DT_histo1->GetZaxis()->SetRangeUser(0,100); 
    DT_histo1->GetXaxis()->SetTitleOffset(x_title_offset); DT_histo1->GetYaxis()->SetTitleOffset(y_title_offset);   DT_histo1->GetZaxis()->SetTitleOffset(z_title_offset);
    DT_histo1->GetXaxis()->SetTitle(xaxislabel.c_str());   DT_histo1->GetYaxis()->SetTitle(yaxislabel.c_str());     DT_histo1->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label1.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "DT");
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c4, pdfname4.c_str());

    y_range_2 = 10.00;
    std::stringstream pdfname5ss; pdfname5ss<<"STAMuonHits_2D_2023_"<<"Percentage"<<"_CSC"; std::string pdfname5 = pdfname5ss.str(); 
    TCanvas * c5 = new TCanvas(pdfname5.c_str(), "", 800, 400); c5->cd(); c5->SetTicks(1,1); c5->SetLeftMargin(0.05); c5->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    CSC_histo1->Draw("textcolz"); /*l7->Draw();*/ CSC_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); CSC_histo1->GetZaxis()->SetRangeUser(0,100);
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, label1.c_str()); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "CSC");
    CSC_histo1->GetXaxis()->SetTitleOffset(x_title_offset); CSC_histo1->GetYaxis()->SetTitleOffset(y_title_offset);   CSC_histo1->GetZaxis()->SetTitleOffset(z_title_offset);
    CSC_histo1->GetXaxis()->SetTitle(xaxislabel.c_str());   CSC_histo1->GetYaxis()->SetTitle(yaxislabel.c_str());     CSC_histo1->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c5, pdfname5.c_str());



}

void savePlot(TCanvas * c, TString name) {
  c->SaveAs(name+".png");
  c->SaveAs(name+".eps");
  gSystem->Exec("epstopdf "+name+".eps");
  std::cout<<"Plot "<<name<<" saved"<<std::endl;
}
