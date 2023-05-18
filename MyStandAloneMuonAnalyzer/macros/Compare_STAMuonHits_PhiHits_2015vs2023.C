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

void Compare_STAMuonHits_PhiHits_2015vs2023() {

    // Title
    // -----
    std::string title = "62Xvs71X";


    // Main Comparison
    // ---------------
    std::string label1 = "CMSSW 62X SLHC (2023 Geometry)"; 
    std::string label2 = "CMSSW 71X (2015 Geometry)"; 

    // Additional Information
    // ----------------------
    std::string info1 = "All #varphi hits of the Stand Alone Muon";
    std::string info2 = "1D DT and 2D CSC hits and Rechits (RPC/GEM)";
    std::string info3 = "Sample used: RelValZmumuJets_Pt_20_300";

    double y_range_1 = 00.00;
    double y_range_2 = 50.00;
    double y_title_offset = 1.20;

    // Reading the filename from command line
    // --------------------------------------
    // std::string RootFileName = gFile->GetName();
    
    // Filename hardcoded
    // ------------------
    std::string RootFileName1 = "STAMuon_62X_2023.root";
    std::string RootFileName2 = "STAMuon_70X_2015.root";

    // Filename base
    // -------------
    std::string Histo_RPCb = "AllHits/Rechits_All_PhiHits_Eta_1D";
    std::string Histo_RPCf = "SegmentsAndHits/Rechits_RPC_Eta_1D";
    std::string Histo_CSC  = "AllHits/Rechits_CSC_Hits_Eta_1D";
    std::string Histo_DT   = "AllHits/Rechits_DT_Hits_PhiSL_Eta_1D";


    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1,0);
    gStyle->SetOptStat(000000);
    gStyle->SetOptTitle(0);
    gStyle->SetGridStyle(1);    // grid
    gROOT->UseCurrentStyle();
    
    TFile * infile1             = new TFile((RootFileName1).c_str(),"READ");     std::cout<<"File 1  = "<<infile1<<std::endl; 
    TFile * infile2             = new TFile((RootFileName2).c_str(),"READ");     std::cout<<"File 2  = "<<infile2<<std::endl; 

    TH1F * RPCb_histo1          = (TH1F*) infile1->Get(Histo_RPCb.c_str());  std::cout<<"Histo  = "<<RPCb_histo1<<std::endl; 
    TH1F * RPCb_histo2          = (TH1F*) infile2->Get(Histo_RPCb.c_str());  std::cout<<"Histo  = "<<RPCb_histo2<<std::endl; 
    TH1F * RPCf_histo1          = (TH1F*) infile1->Get(Histo_RPCf.c_str());  std::cout<<"Histo  = "<<RPCf_histo1<<std::endl; 
    TH1F * RPCf_histo2          = (TH1F*) infile2->Get(Histo_RPCf.c_str());  std::cout<<"Histo  = "<<RPCf_histo2<<std::endl; 
    TH1F * CSC_histo1           = (TH1F*) infile1->Get(Histo_CSC.c_str());   std::cout<<"Histo  = "<<CSC_histo1<<std::endl; 
    TH1F * CSC_histo2           = (TH1F*) infile2->Get(Histo_CSC.c_str());   std::cout<<"Histo  = "<<CSC_histo2<<std::endl; 
    TH1F * DT_histo1            = (TH1F*) infile1->Get(Histo_DT.c_str());    std::cout<<"Histo  = "<<DT_histo1<<std::endl; 
    TH1F * DT_histo2            = (TH1F*) infile2->Get(Histo_DT.c_str());    std::cout<<"Histo  = "<<DT_histo2<<std::endl; 

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
    std::string Yaxis = "Average \# Hits";

    RPCb_histo1->SetTitle("All");  RPCb_histo1->GetXaxis()->SetTitle(Xaxis.c_str()); RPCb_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    RPCf_histo1->SetTitle("RPC");  RPCf_histo1->GetXaxis()->SetTitle(Xaxis.c_str()); RPCf_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    CSC_histo1->SetTitle("CSC");   CSC_histo1->GetXaxis()->SetTitle(Xaxis.c_str());  CSC_histo1->GetYaxis()->SetTitle(Yaxis.c_str());
    DT_histo1->SetTitle("DT");     DT_histo1->GetXaxis()->SetTitle(Xaxis.c_str());   DT_histo1->GetYaxis()->SetTitle(Yaxis.c_str());

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

    std::stringstream pdfname1ss; pdfname1ss<<"Compare_STAMuonHits_"<<title<<"_AllPlots"; std::string pdfname1 = pdfname1ss.str(); 
    TCanvas * c1 = new TCanvas(pdfname1.c_str(), "", 800, 600);
    
    c1->Divide(2,2);
    c1->cd(1)->SetTicks(1,1);
    c1->cd(2)->SetTicks(1,1);
    c1->cd(3)->SetTicks(1,1);
    c1->cd(4)->SetTicks(1,1);

    double infoX_x_right = 0.5; 
    double info1_y = 0.80, info2_y = 0.75, info3_y = 0.70;
    y_range_2 = 60.00;
    c1->cd(1); RPCb_histo1->Draw("HF"); RPCb_histo2->Draw("HFsame"); l1->Draw(); RPCb_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); RPCb_histo1->GetYaxis()->SetTitleOffset(y_title_offset);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 8 TeV"); latex_cmslab.DrawLatex(0.85, 0.925,"All"); gPad->RedrawAxis();
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    y_range_2 = 10.00;
    c1->cd(2); RPCf_histo1->Draw("HF"); RPCf_histo2->Draw("Hsame"); l2->Draw(); RPCf_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); RPCf_histo1->GetYaxis()->SetTitleOffset(y_title_offset);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 8 TeV"); latex_cmslab.DrawLatex(0.85, 0.925,"RPC"); gPad->RedrawAxis();
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    y_range_2 = 50.00;
    c1->cd(3); DT_histo1->Draw("HF");   DT_histo2->Draw("HFsame");   l4->Draw(); DT_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); DT_histo1->GetYaxis()->SetTitleOffset(y_title_offset);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 8 TeV"); latex_cmslab.DrawLatex(0.85, 0.925,"DT"); gPad->RedrawAxis();
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    y_range_2 = 50.00;
    c1->cd(4); CSC_histo1->Draw("HF");  CSC_histo2->Draw("HFsame");  l3->Draw(); CSC_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); CSC_histo1->GetYaxis()->SetTitleOffset(y_title_offset);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 8 TeV"); latex_cmslab.DrawLatex(0.85, 0.925,"CSC"); gPad->RedrawAxis();
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    c1->Modified();
    c1->Update();

    // Save this Canvas
    savePlot(c1, pdfname1.c_str());


    TLegend * l5 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l5->SetLineColor(1);    l5->SetLineStyle(1);  l5->SetLineWidth(1);  l5->SetFillColor(4000); l5->SetBorderSize(1);
    TLegend * l6 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l6->SetLineColor(1);    l6->SetLineStyle(1);  l6->SetLineWidth(1);  l6->SetFillColor(4000); l6->SetBorderSize(1);
    TLegend * l7 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l7->SetLineColor(1);    l7->SetLineStyle(1);  l7->SetLineWidth(1);  l7->SetFillColor(4000); l7->SetBorderSize(1);
    TLegend * l8 = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, NULL,"brNDC"); l8->SetLineColor(1);    l8->SetLineStyle(1);  l8->SetLineWidth(1);  l8->SetFillColor(4000); l8->SetBorderSize(1);

    l5->AddEntry(RPCb_histo1, RPC_b_1_label.c_str(), "F");     l5->AddEntry(RPCb_histo2, RPC_b_2_label.c_str(), "F");
    l6->AddEntry(RPCf_histo1, RPC_f_1_label.c_str(), "F");     l6->AddEntry(RPCf_histo2, RPC_f_2_label.c_str(), "F");
    l7->AddEntry(CSC_histo1,  CSC_1_label.c_str(),   "F");     l7->AddEntry(CSC_histo2,  CSC_2_label.c_str(),   "F");
    l8->AddEntry(DT_histo1,   DT_1_label.c_str(),    "F");     l8->AddEntry(DT_histo2,   DT_2_label.c_str(),    "F");

    // Save the Other Canvasses
    std::stringstream pdfname2ss; pdfname2ss<<"Compare_STAMuonHits_PhiHits_"<<title<<"_All"; std::string pdfname2 = pdfname2ss.str(); 
    TCanvas * c2 = new TCanvas(pdfname2.c_str(), "", 800, 600); c2->cd(); c2->SetTicks(1,1);
    y_range_2 = 60.00;
    RPCb_histo1->Draw("HF"); RPCb_histo2->Draw("HFsame"); l5->Draw(); RPCb_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); RPCb_histo1->GetYaxis()->SetTitleOffset(y_title_offset);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.85, 0.925,"All"); gPad->RedrawAxis();
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c2, pdfname2.c_str());

    std::stringstream pdfname3ss; pdfname3ss<<"Compare_STAMuonHits_PhiHits_"<<title<<"_RPC"; std::string pdfname3 = pdfname3ss.str(); 
    TCanvas * c3 = new TCanvas(pdfname3.c_str(), "", 800, 600); c3->cd(); c3->SetTicks(1,1);
    y_range_2 = 10.00;
    RPCf_histo1->Draw("HF"); RPCf_histo2->Draw("Hsame"); l6->Draw(); RPCf_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); RPCf_histo1->GetYaxis()->SetTitleOffset(y_title_offset);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.85, 0.925,"RPC"); gPad->RedrawAxis();
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c3, pdfname3.c_str());

    std::stringstream pdfname4ss; pdfname4ss<<"Compare_STAMuonHits_PhiHits_"<<title<<"_DT"; std::string pdfname4 = pdfname4ss.str(); 
    TCanvas * c4 = new TCanvas(pdfname4.c_str(), "", 800, 600); c4->cd(); c4->SetTicks(1,1);
    y_range_2 = 50.00;
    DT_histo1->Draw("HF");   DT_histo2->Draw("HFsame");   l8->Draw(); DT_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); DT_histo1->GetYaxis()->SetTitleOffset(y_title_offset);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.85, 0.925,"DT"); gPad->RedrawAxis();
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c4, pdfname4.c_str());

    std::stringstream pdfname5ss; pdfname5ss<<"Compare_STAMuonHits_PhiHits_"<<title<<"_CSC"; std::string pdfname5 = pdfname5ss.str(); 
    TCanvas * c5 = new TCanvas(pdfname5.c_str(), "", 800, 600); c5->cd(); c5->SetTicks(1,1);
    y_range_2 = 50.00;
    CSC_histo1->Draw("HF");  CSC_histo2->Draw("HFsame");  l7->Draw(); CSC_histo1->GetYaxis()->SetRangeUser(y_range_1,y_range_2); CSC_histo1->GetYaxis()->SetTitleOffset(y_title_offset);
    latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 13 TeV"); latex_cmslab.DrawLatex(0.85, 0.925,"CSC"); gPad->RedrawAxis();
    latex_info.DrawLatex(infoX_x_right, info1_y, info1.c_str()); latex_info.DrawLatex(infoX_x_right, info2_y, info2.c_str()); latex_info.DrawLatex(infoX_x_right, info3_y, info3.c_str());
    savePlot(c5, pdfname5.c_str());



}

void savePlot(TCanvas * c, TString name) {
  c->SaveAs(name+".png");
  c->SaveAs(name+".eps");
  gSystem->Exec("epstopdf "+name+".eps");
  std::cout<<"Plot "<<name<<" saved"<<std::endl;
}
