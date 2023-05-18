#ifndef SGraphContainer_h
#define SGraphContainer_h

#include "string"
#include "map"
#include "iostream"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TDirectoryFile.h"

class SGraphContainer
{
 public:
  SGraphContainer();
  ~SGraphContainer();
  void book(std::string path,TGraph*h);
  void book(TGraph*);//same as book("",h)

  TGraph * get(std::string path,std::string name);
  TGraph * get(std::string name);//same as get("",name)
  TGraphErrors * getE(std::string path,std::string name);
  TGraphErrors * getE(std::string name);//same as get("",name) 
  TGraphAsymmErrors * getA(std::string path,std::string name);
  TGraphAsymmErrors * getA(std::string name);//same as get("",name) 

  void write(TDirectoryFile * d);
 protected:
  std::map<std::string,std::map<std::string,TGraph*> > graph_;
};

SGraphContainer::SGraphContainer(){}
SGraphContainer::~SGraphContainer(){}

void SGraphContainer::book(std::string path, TGraph * h) {
  std::map<std::string,std::map<std::string, TGraph *> >::iterator dirit = graph_.find(path);
  if(dirit == graph_.end()) {
      std::pair<std::string,std::map<std::string, TGraph *> > thepair(path,std::map<std::string,TGraph*>());
      dirit = graph_.insert(thepair).first;
  }
  std::map<std::string, TGraph *>::iterator graphit = dirit->second.insert(std::pair<std::string,TGraph*>(h->GetName(),h)).first;
}

void SGraphContainer::book(TGraph*h) { this->book("",h); }

TGraphAsymmErrors * SGraphContainer::getA(std::string path,std::string name) { return (TGraphAsymmErrors*)this->get(path,name); }
TGraphErrors      * SGraphContainer::getE(std::string path,std::string name) { return (TGraphErrors*)     this->get(path,name); }
TGraph            * SGraphContainer::get( std::string path,std::string name) {
  std::map<std::string,std::map<std::string,TGraph *> >::iterator dirit = graph_.find(path);
  if(dirit == graph_.end())          return 0; 
  std::map<std::string,TGraph *>::iterator graphit = dirit->second.find(name);
  if(graphit == dirit->second.end()) return 0;
  else                               return graphit->second;
}






TGraphAsymmErrors * SGraphContainer::getA(std::string name) { return (TGraphAsymmErrors*)this->get("",name); }
TGraphErrors      * SGraphContainer::getE(std::string name) { return (TGraphErrors*)     this->get("",name); }
TGraph            * SGraphContainer::get(std::string name)  { return                      this->get("",name); }

void SGraphContainer::write(TDirectoryFile * d) {
  for(std::map<std::string,std::map<std::string,TGraph*> >::iterator dirit = graph_.begin();dirit!=graph_.end();++dirit) {
    TDirectoryFile * curdir = d;
    if(dirit->first!="") {
      std::vector<std::string> path;
      path.reserve(5);
      std::string str=dirit->first;
      unsigned int pos = str.find_first_of('/');
      while(pos != std::string::npos) {
	std::string element = str.substr(0,pos);
	path.push_back(element);
	str = str.substr(pos+1);
	pos = str.find_first_of('/');
      }
      if(str!="")
	path.push_back(str);
      for(unsigned int i=0; i<path.size(); ++i) {
	TDirectoryFile * curdir_ = (TDirectoryFile*)curdir->Get(path[i].c_str());
	if(!curdir_) curdir = (TDirectoryFile*)curdir->mkdir(path[i].c_str());
	else         curdir = curdir_;
      }
    }
    curdir->cd();
    for(std::map<std::string,TGraph*>::iterator graphit = dirit->second.begin();graphit!=dirit->second.end();++graphit) {
      std::string graphname = graphit->second->GetName();
      graphit->second->Write(graphname.c_str());
    }
  }
}

#endif
