#ifndef SHistContainer_h
#define SHistContainer_h

#include "string"
#include "map"
#include "iostream"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TDirectoryFile.h"

class SHistContainer
{
 public:
  SHistContainer();
  ~SHistContainer();
  void book(std::string path,TH1*h);
  void book(TH1*);//same as book("",h)
  TH1 * get(std::string path,std::string name);
  TH1 * get(std::string name);//same as get("",name)
  TH2 * get2(std::string path,std::string name);
  TH2 * get2(std::string name);//same as get("",name)
  TProfile * getP(std::string path,std::string name);
  TProfile * getP(std::string name);//same as get("",name)
  void write(TDirectoryFile * d);
 protected:
  std::map<std::string,std::map<std::string,TH1*> > hist_;
};

SHistContainer::SHistContainer(){}
SHistContainer::~SHistContainer(){}

void SHistContainer::book(std::string path, TH1 * h) {
  std::map<std::string,std::map<std::string, TH1 *> >::iterator dirit = hist_.find(path);
  if(dirit == hist_.end()) {
      std::pair<std::string,std::map<std::string, TH1 *> > thepair(path,std::map<std::string,TH1*>());
      dirit = hist_.insert(thepair).first;
  }
  std::map<std::string, TH1 *>::iterator histit = dirit->second.insert(std::pair<std::string,TH1*>(h->GetName(),h)).first;
}

void SHistContainer::book(TH1*h) { this->book("",h); }

TProfile * SHistContainer::getP(std::string path,std::string name) { return (TProfile*)this->get(path,name); }

TH2 * SHistContainer::get2(std::string path,std::string name) { return (TH2*)this->get(path,name); }

TH1 * SHistContainer::get(std::string path,std::string name) {
  std::map<std::string,std::map<std::string,TH1 *> >::iterator dirit = hist_.find(path);
  if(dirit == hist_.end())          return 0;
  std::map<std::string,TH1 *>::iterator histit = dirit->second.find(name);
  if(histit == dirit->second.end()) return 0;
  else                              return histit->second;
}

TProfile * SHistContainer::getP(std::string name) { return (TProfile*)this->get("",name); }

TH2 * SHistContainer::get2(std::string name) { return (TH2*)this->get("",name); }

TH1 * SHistContainer::get(std::string name)  { return this->get("",name); }

void SHistContainer::write(TDirectoryFile * d) {
  for(std::map<std::string,std::map<std::string,TH1*> >::iterator dirit = hist_.begin();dirit!=hist_.end();++dirit) {
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
    for(std::map<std::string,TH1*>::iterator histit = dirit->second.begin();histit!=dirit->second.end();++histit)
      histit->second->Write();
  }
}

#endif
