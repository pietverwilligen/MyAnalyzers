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

void Investigate_PT_StaMuonHits_2015() {

    // Title
    // -----
    std::string title = "Stand Alone Muons Investigation";


    // Main Comparison
    // ---------------
    std::string label1 = "All Stand Alone Muons (Sta Mu)"; 
    std::string label2 = "Sta Mu (1.7 < |#eta| < 2.5 && 2 const hits)"; 

    // Additional Information
    // ----------------------
    std::string info1 = "Hits used in Stand Alone Muon Fit";
    std::string info2 = "these are Segments (CSC,DT) and Rechits (RPC)";
    std::string info3 = "Sample used: RelValZmumuJets_Pt_20_300";

    // double y_range_1 = 00.00;
    // double y_range_2 = 01.00;
    double y_range_1 = 1E-5;
    double y_range_2 = 1E0;
    double x_title_offset = 0.75;
    double y_title_offset = 0.75;
    double z_title_offset = 0.75;

    std::string xaxislabel = "STA Muon p_{T}";
    std::string yaxislabel = "Fraction";
    // std::string zaxislabel = "Percentage";

    // Reading the filename from command line
    // --------------------------------------
    // std::string RootFileName = gFile->GetName();
    
    // Filename hardcoded
    // ------------------
    std::string RootFileName1 = "STAMuon_62X_2023.root";
    std::string RootFileName2 = "STAMuon_70X_2015.root";

    // Filename base
    // -------------
    std::string Histo_StandAlone_PT_All = "StandAloneMuon_PT";
    std::string Histo_StandAlone_PT_Sel = "StandAloneMuon_1p8To2p5_2Hit_PT";


    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1,0);
    gStyle->SetOptStat(000000);
    gStyle->SetOptTitle(0);
    gStyle->SetGridStyle(1);    // grid
    gROOT->UseCurrentStyle();

    // TFile * infile1             = new TFile((RootFileName1).c_str(),"READ");     std::cout<<"File 1  = "<<infile1<<std::endl;     
    TFile * infile2             = new TFile((RootFileName2).c_str(),"READ");     std::cout<<"File 2  = "<<infile2<<std::endl; 


    TH2F * Sta_histo1          = (TH2F*) infile2->Get(Histo_StandAlone_PT_All.c_str());  std::cout<<"Histo  = "<<Histo_StandAlone_PT_All<<std::endl; 
    TH2F * Sta_histo2          = (TH2F*) infile2->Get(Histo_StandAlone_PT_Sel.c_str());  std::cout<<"Histo  = "<<Histo_StandAlone_PT_Sel<<std::endl; 


    Sta_histo1->SetLineColor(kBlue);  Sta_histo1->SetFillColor(kBlue);  Sta_histo1->SetFillStyle(3005);
    Sta_histo1->SetLineWidth(1);

    Sta_histo2->SetLineColor(kRed);   Sta_histo2->SetFillColor(kRed); Sta_histo2->SetFillStyle(3004);
    Sta_histo2->SetLineWidth(1);

    std::string Xaxis = "p_{T}";
    std::string Yaxis = "Fraction";

    Sta_histo1->SetTitle("All");  Sta_histo1->GetXaxis()->SetTitle(Xaxis.c_str()); Sta_histo1->GetYaxis()->SetTitle(Yaxis.c_str());

    int Sta_1_entries = Sta_histo1->GetEntries();                    int Sta_2_entries = Sta_histo2->GetEntries();
    std::stringstream ss_Sta_1_label; ss_Sta_1_label << label1;      std::stringstream ss_Sta_2_label; ss_Sta_2_label << label2;
    std::string Sta_1_label = ss_Sta_1_label.str();                  std::string Sta_2_label = ss_Sta_2_label.str();

    Sta_histo1->Scale(1.0/Sta_1_entries);     Sta_histo2->Scale(1.0/Sta_2_entries); 


    double leg_x1 = 0.55, leg_x2 = 0.85, leg_y1 = 0.700, leg_y2 = 0.850; 
    TLegend * l1 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l1->SetLineColor(1);    l1->SetLineStyle(1);  l1->SetLineWidth(1);  l1->SetFillColor(4000); l1->SetBorderSize(1);
    l1->AddEntry(Sta_histo1, Sta_1_label.c_str(), "F");     l1->AddEntry(Sta_histo2, Sta_2_label.c_str(), "F");

    double leg_x3 = 0.15, leg_x4 = 0.45, leg_y3 = 0.150, leg_y4 = 0.300; 
    TLegend * l2 = new TLegend(leg_x3, leg_y3, leg_x4, leg_y4, NULL,"brNDC"); l2->SetLineColor(1);    l2->SetLineStyle(1);  l2->SetLineWidth(1);  l2->SetFillColor(4000); l2->SetBorderSize(1);
    l2->AddEntry(Sta_histo1, Sta_1_label.c_str(), "F");     l2->AddEntry(Sta_histo2, Sta_2_label.c_str(), "F");

    TLatex latex_cmslab; latex_cmslab.SetNDC(); latex_cmslab.SetTextSize(0.035); latex_cmslab.SetTextAlign(11);
    TLatex latex_info;   latex_info.SetNDC();   latex_info.SetTextSize(0.025);    latex_info.SetTextAlign(31);

    // double infoX_x_right = 0.6275; 
    // double info1_y = 0.80, info2_y = 0.75, info3_y = 0.70;

    y_range_2 = 0.30;    
    std::stringstream pdfname1ss; pdfname1ss<<"STAMuonPtSpectrum_2015_"<<"StaMuInvestigation"<<"_Lin"; std::string pdfname1 = pdfname1ss.str(); 
    TCanvas * c1 = new TCanvas(pdfname1.c_str(), "", 600, 400); c1->cd(); c1->SetTicks(1,1); c1->SetLeftMargin(0.05); c1->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    Sta_histo1->Draw(); Sta_histo2->Draw("same");           Sta_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); // Sta_histo2->GetZaxis()->SetRangeUser(0,100); 
    Sta_histo1->GetXaxis()->SetTitleOffset(x_title_offset); Sta_histo1->GetYaxis()->SetTitleOffset(y_title_offset);    // Sta_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    Sta_histo1->GetXaxis()->SetTitle(xaxislabel.c_str());   Sta_histo1->GetYaxis()->SetTitle(yaxislabel.c_str());      // Sta_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, "Stand Alone Muons"); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "All");
    // latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    l1->Draw();
    c1->Update();
    savePlot(c1, pdfname1.c_str());

    y_range_2 = 1.00;
    std::stringstream pdfname2ss; pdfname2ss<<"STAMuonPtSpectrum_2015_"<<"StaMuInvestigation"<<"_Log"; std::string pdfname2 = pdfname2ss.str(); 
    TCanvas * c2 = new TCanvas(pdfname2.c_str(), "", 600, 400); c2->cd(); c2->SetTicks(1,1); c2->SetLeftMargin(0.05); c2->SetRightMargin(0.10); gStyle->SetPaintTextFormat("3.0f");
    Sta_histo1->Draw(); Sta_histo2->Draw("same");           Sta_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); // Sta_histo2->GetZaxis()->SetRangeUser(0,100); 
    Sta_histo1->GetXaxis()->SetTitleOffset(x_title_offset); Sta_histo1->GetYaxis()->SetTitleOffset(y_title_offset);    // Sta_histo2->GetZaxis()->SetTitleOffset(z_title_offset);
    Sta_histo1->GetXaxis()->SetTitle(xaxislabel.c_str());   Sta_histo1->GetYaxis()->SetTitle(yaxislabel.c_str());      // Sta_histo2->GetZaxis()->SetTitle(zaxislabel.c_str());
    latex_cmslab.DrawLatex(0.05, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.70, 0.925, "Stand Alone Muons"); gPad->RedrawAxis(); latex_cmslab.DrawLatex(0.125, 0.85, "All");
    // latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    l2->Draw();
    c2->SetLogy();
    c2->SetLogx();
    c2->Update();
    savePlot(c2, pdfname2.c_str());


}

void savePlot(TCanvas * c, TString name) {
  c->SaveAs(name+".png");
  c->SaveAs(name+".eps");
  gSystem->Exec("epstopdf "+name+".eps");
  std::cout<<"Plot "<<name<<" saved"<<std::endl;
}
