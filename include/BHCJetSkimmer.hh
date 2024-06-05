#ifndef BHCJETSKIMMER_HH
#define BHCJETSKIMMER_HH
#include "JetSimProducer.hh"
#include "BaseSkimmer.hh"
#include "BaseTree.hh"
#include "TGraph.h"

using node = BaseTree::node;
using procCat = BaseSkimmer::procCat;

class BHCJetSkimmer{
	public:
		BHCJetSkimmer(){
			_evti = 0;
			_evtj = 0;
			_gev = 1./10.;
			_oname = "";
			_radius = 0;
			_smear = true;
			_timesmear = false;
		}

		virtual ~BHCJetSkimmer(){ }

		BHCJetSkimmer(TFile* file){
			_prod = new JetSimProducer(file);

			_oname = "";
			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_gev = 1./10.;
			_radius = 0;
			_smear = true;
			_timesmear = false;
				
	
			graphs.push_back(nrhs_comptime);

			_hists1D.push_back(nClusters);
			_hists1D.push_back(nSubClusters);	
			_hists1D.push_back(predJet_subClusterEnergy);
			_hists1D.push_back(predJet_subClusterEtaCenter);
			_hists1D.push_back(predJet_subClusterPhiCenter);
			_hists1D.push_back(predJet_subClusterTimeCenter);
			_hists1D.push_back(predJet_dR);
			_hists1D.push_back(predJet_energy);
			_hists1D.push_back(predJet_pt);
			_hists1D.push_back(predJet_mass);
			_hists1D.push_back(jetGenE_sigmaDeltaPt_predGen);
			_hists1D.push_back(predGen_nJets);
			_hists1D.push_back(predJet_subClusterEtaSigma);
			_hists1D.push_back(predJet_subClusterPhiSigma);
			_hists1D.push_back(predJet_subClusterTimeSigma);
			_hists1D.push_back(predJet_subClusteretaPhiCov);
			_hists1D.push_back(predJet_subClustertimeEtaCov);
			_hists1D.push_back(predJet_subClustertimePhiCov);
			_hists1D.push_back(nRecoJets);
			_hists1D.push_back(recoJet_dR);
			_hists1D.push_back(recoJet_energy);
			_hists1D.push_back(recoJet_pt);
			_hists1D.push_back(recoJet_mass);
			_hists1D.push_back(jetGenE_sigmaDeltaPt_recoGen);
			_hists1D.push_back(recoGen_nJets);
			_hists1D.push_back(recoGen_jetPtRatio);

			_hists2D.push_back(jetGenE_diffDeltaPt_predGen);
			_hists2D.push_back(jetGenE_diffDeltaPt_recoGen);
			_hists2D.push_back(genPt_recoPt);

		}
		void Skim();
		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else return; 
		}
	

		void CleanTrees(const vector<node*>& trees){
			_trees.clear();
			for(int i = 0; i < trees.size(); i++){
				if(trees[i] == nullptr) continue;
				//check for mirrored point - would be double counted
				if(trees[i]->points->mean().at(1) > 2*acos(-1) || trees[i]->points->mean().at(1) < 0) continue;	
				_trees.push_back(trees[i]);
			}
		}


		void TreesToJets(){
			_predJets.clear();
			vector<JetPoint> rhs;
			double x, y, z, eta, phi, t, theta;
			BayesPoint vertex({_pvx, _pvy, _pvz});

			for(int i = 0; i < _trees.size(); i++){
				//get points from tree
				PointCollection* pc = _trees[i]->points;
				rhs.clear();
				//loop over points
				//pc->Print();
				double jeta = 0;
				double jphi = 0;
				double je = 0;
				for(int p = 0; p < pc->GetNPoints(); p++){
					eta = pc->at(p).at(0);
					phi = pc->at(p).at(1);
					t = pc->at(p).at(2);
					//translate eta, phi to x, y, z
					theta = 2*atan2(1,exp(eta));
					x = _radius*cos(phi);
					y = _radius*sin(phi);
					z = _radius/tan(theta);	
					
					//for checking calculations
					double rtheta = atan2(_radius,z);
					//if(phoz < 0) rtheta = atan2(129.,phoz)+acos(-1)/2.;
					double reta = -log(tan(rtheta/2));
					double rphi = atan2(y,x);
		
					//cout << "eta: " << eta << " reta: " << reta << " phi: " << phi << " rphi: " << rphi << " theta: " << theta << " rtheta: " << rtheta << endl;
					
					//declare JetPoint with x, y, z, t
					JetPoint jp(x, y, z, t);
					//cout << "jp eta " << jp.eta() << " eta " << eta << " jp phi " << jp.phi() << " phi " << phi << endl;
					//add JetPoint to list of rhs
					jp.SetEnergy(pc->at(p).w()/_gev);
					jp.SetWeight(pc->at(p).w());
					rhs.push_back(jp);
					jeta += eta;//*pc->at(p).w()/_gev;
					jphi += phi;//*pc->at(p).w()/_gev;
					je += pc->at(p).w()/_gev;
				}
				//create new Jet
				Jet predJet(rhs, BayesPoint({_pvx, _pvy, _pvz}));
				//set PV info
				//predJet.SetVertex(Point({_pvx,_pvy,_pvz}));
				//set constituents (subclusters) here with model from tree
				int nsubclusters = _trees[i]->model->GetNClusters();
				double Ek; //effective energy of constituent
				vector<double> norms;
				_trees[i]->model->GetNorms(norms);
				map<string,Matrix> params;
				for(int k = 0; k < nsubclusters; k++){
					params = _trees[i]->model->GetPriorParameters(k);
					Ek = norms[k]/_gev;
					predJet.AddConstituent(params,Ek);
				}
				//TODO: set covariance and mean for jets with 1+ subcluster
				if(nsubclusters < 2){
					params = _trees[i]->model->GetPriorParameters(0);
					predJet.SetCovariance(params["cov"]);
					predJet.SetCenter(params["mean"]);
				}
				//cout << "jet mean eta " << jeta/je << " mean phi " << jphi/je << " total e " << je << endl;
				//cout << "jet mean eta " << jeta/(double)pc->GetNPoints() << " mean phi " << jphi/(double)pc->GetNPoints() << " total e " << je << " " << predJet.E() << " mass " << predJet.mass() << endl;
				//pc->Print();
				//add Jet to jets	
				_predJets.push_back(predJet);	
				//vector<Jet> babies = _predJets[_predJets.size()-1].GetConstituents();
				//vector<JetPoint> rrhs = _predJets[_predJets.size()-1].GetJetPoints();
			}

		}


		/*
		void MatchJetsToTracks(){
			int nTracks = _base->Track_px->size();
			cout << "nTracks " << nTracks << endl;
			int nConsts; //# of constituents
			for(int j = 0; j < _predJets.size(); j++){
				vector<Jet>& consts = _predJets[j].GetConstituents();
				for(int c = 0; c < consts.size(); c++){
					//consider putting a minimum dR cut on subcluster/track matching
					double drmin = 999;
					double dr;
					int bestIdx;
					for(int t = 0; t < nTracks; t++){
						//find min dR bw constituent and track
						//may want to do a momentum or energy matching too
						dr = dR(consts[c].eta(), consts[c].phi(), _base->Track_eta->at(t), _base->Track_phi->at(t));
						if(dr < drmin){ drmin = dr; bestIdx = t; }
					}
					cout << "jet #" << j << " eta " << _predJets[j].eta() << " phi " << _predJets[j].phi() << " E " << _predJets[j].E() << " constituent #" << c << " eta " << consts[c].eta() << " phi " << consts[c].phi() << " E " << consts[c].E() << " best track match track #" << bestIdx << " with dR " << drmin << endl;
					//cout << "jet #" << j << " eta " << _predJets[j].eta() << " phi " << _predJets[j].phi() << " E " << _predJets[j].E() << " constituent #" << c << " eta " << subcl.eta() << " phi " << subcl.phi() << " E " << subcl.E() << " best track match track #" << bestIdx << " with dR " << drmin << endl;
					//set momentum of constituent based on track momentum
					consts[c].SetP(_base->Track_px->at(bestIdx), _base->Track_py->at(bestIdx), _base->Track_pz->at(bestIdx));
					//cout << "consts #" << c << " px " << consts[c].px() << " py " << consts[c].py() << " pz " << consts[c].pz() << endl;
				}
				//set momentum of jet to sum of constituents
				cout << "Nconsts " << _predJets[j].GetNConstituents() << endl;
				_predJets[j].SetP();
				cout << "jet #" << j << " px " << _predJets[j].px() << " py " << _predJets[j].py() << " pz " << _predJets[j].pz() << " mass " << _predJets[j].mass() << " energy " << _predJets[j].E() << endl;
			}


		}
		*/
	
		void FillPredJetHists(){
			int njets;
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				njets = _predJets.size();
				for(int j = 0; j < _predJets.size(); j++){
					//cout << "jet #" << j << " phi " << _predJets[j].phi() << " eta " << _predJets[j].eta() << " energy " << _predJets[j].E() <<  " mass " << _predJets[j].mass() << " nConstituents " << _predJets[j].GetNConstituents() << " nRhs " << _predJets[j].GetNRecHits() << " pt " << _predJets[j].pt() << endl;
					_procCats[p].hists1D[0][7]->Fill(_predJets[j].e());
					_procCats[p].hists1D[0][8]->Fill(_predJets[j].pt());
					_procCats[p].hists1D[0][9]->Fill(_predJets[j].mass());
					//fill dR hist - max dR is max dR bw rhs or sub clusters?
					//if subclusters - essentially matching subclusters to particles
					//if rhs - no notion of particles	
					//use subclusters for now (can change to rhs later)
					_procCats[p].hists1D[0][6]->Fill(CalcDr(_predJets[j]));	
				}
				_procCats[p].hists1D[0][0]->Fill(njets);
				//cout << "# pred jets - # gen jets " << njets - (int)_genjets.size() << endl;
				_procCats[p].hists1D[0][11]->Fill(njets - (int)_genjets.size());
			}
		}
		
		void FillRecoJetHists(){
			int njets;
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				njets = _recojets.size();
				for(int j = 0; j < _recojets.size(); j++){
					//cout << "jet #" << j << " phi " << _recojets[j].phi() << " eta " << _recojets[j].eta() << " energy " << _recojets[j].E() <<  " mass " << _recojets[j].mass() << " nConstituents " << _recojets[j].GetNConstituents() << " nRhs " << _recojets[j].GetNRecHits() << " pt " << _recojets[j].pt() << endl;
					_procCats[p].hists1D[0][20]->Fill(_recojets[j].e());
					//dr hist is #19
					_procCats[p].hists1D[0][21]->Fill(_recojets[j].pt());
					_procCats[p].hists1D[0][22]->Fill(_recojets[j].mass());
				}
				_procCats[p].hists1D[0][18]->Fill(njets);
				//cout << "hist name " << _procCats[p].hists1D[0][18]->GetName() << " nentries " << _procCats[p].hists1D[0][18]->GetEntries() << endl;
				//cout << "# pred jets - # gen jets " << njets - (int)_genjets.size() << endl;
				_procCats[p].hists1D[0][24]->Fill(njets - (int)_genjets.size());
			}
		}

		void FillResolutionHists(){
			//need to gen match jets to find difference in pt bw reco - gen
			//dr match (maybe do dE match later)
			for(int p = 0; p < _procCats.size(); p++){
				double dr;
				int bestIdx = 0;
				//cout << "proc " << p << " # reco jets " << _recojets.size() << " # gen jets " << _genjets.size() << endl;
				//reco jets
				for(int j = 0; j < _recojets.size(); j++){
					dr = 999;
					for(int g = 0; g < _genjets.size(); g++){
						if(dR(_recojets[j].eta(), _recojets[j].phi(), _genjets[g].eta(), _genjets[g].phi()) < dr){
							dr = dR(_recojets[j].eta(), _recojets[j].phi(), _genjets[g].eta(), _genjets[g].phi());
							bestIdx = g;
						}
					}
					//cout << "jet #" << j << " has best match with gen jet #" << bestIdx << " with dr " << dr << " reco E " << _predJets[j].E() << " gen energy " << _genjets[bestIdx].E() << " reco pt " << _recojets[j].pt() << " gen pt " << _genjets[bestIdx].pt() << endl;
					_procCats[p].hists2D[0][1]->Fill(_genjets[bestIdx].E(), _recojets[j].pt() - _genjets[bestIdx].pt());
					_procCats[p].hists1D[0][25]->Fill(_recojets[j].pt()/_genjets[bestIdx].pt());
					_procCats[p].hists2D[0][2]->Fill(_genjets[bestIdx].pt(), _recojets[j].pt());
				}
				//predicted jets
				for(int j = 0; j < _predJets.size(); j++){
					dr = 999;
					for(int g = 0; g < _genjets.size(); g++){
						if(dR(_predJets[j].eta(), _predJets[j].phi(), _genjets[g].eta(), _genjets[g].phi()) < dr){
							dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genjets[g].eta(), _genjets[g].phi());
							bestIdx = g;
						}
					}
					//cout << "jet #" << j << " has best match with gen jet #" << bestIdx << " with dr " << dr << " reco E " << _predJets[j].E() << " gen energy " << _genjets[bestIdx].E() << endl;
					_procCats[p].hists2D[0][0]->Fill(_genjets[bestIdx].E(), _predJets[j].pt() - _genjets[bestIdx].pt());
				}
			}


		}

		void FillModelHists(){
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				for(int i = 0; i < _trees.size(); i++){
					//keep checks in just in case
					if(_trees[i] == nullptr) continue;
					//check for mirrored point - would be double counted
					if(_trees[i]->points->mean().at(1) > 2*acos(-1) || _trees[i]->points->mean().at(1) < 0) continue;
					FillModelHists(_trees[i]->model, p);
				}
			}

		}
		
		//all hists referenced here are in hists1D
		void FillModelHists(BasePDFMixture* model, int p){
			map<string, Matrix> params;
			vector<double> eigenvals, norms;
			vector<Matrix> eigenvecs;
			double theta, phi, r, id, npts, E_k;

			int nclusters = model->GetNClusters();
			
			_procCats[p].hists1D[0][1]->Fill(nclusters);
			model->GetNorms(norms);
	
			//center is mm weighted avg of subclusters
			double ceta = 0;
			double cphi = 0;
			double ctime = 0;
			double norm = 0;
			double Etot = 0;
			PointCollection* pts = model->GetData();

	
			//k clusters = k jets in event -> subclusters are mixture model components
			for(int k = 0; k < nclusters; k++){
				E_k = norms[k]/_gev;
				Etot += E_k;

				params = model->GetPriorParameters(k);
				ceta += params["pi"].at(0,0)*params["mean"].at(0,0);
				cphi += params["pi"].at(0,0)*params["mean"].at(1,0);
				//cout << "pi " << params["pi"].at(0,0) << " ctime " << params["mean"].at(0,2) << " cphi " << params["mean"].at(0,1) << " ceta " << params["mean"].at(0,0) << endl;
				//cout << "switched indices - pi " << params["pi"].at(0,0) << " ctime " << params["mean"].at(2,0) << " cphi " << params["mean"].at(1,0) << " ceta " << params["mean"].at(0,0) << endl;
				ctime += params["pi"].at(0,0)*params["mean"].at(2,0);
				norm += params["pi"].at(0,0);
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//total cluster energy
			}
			_procCats[p].hists1D[0][2]->Fill(Etot);
			_procCats[p].hists1D[0][3]->Fill(ceta/norm);
			_procCats[p].hists1D[0][4]->Fill(cphi/norm);
			_procCats[p].hists1D[0][5]->Fill(ctime/norm);
		}
	
		//use subclusters
		//can change to rhs later
		double CalcDr(const Jet& jet){
			double dr = 0;
			double maxDr = 0;
			//if 1 subcluster (should be >1 but idk could happen ig)
			//dR =  sqrt(sigeta*sigeta + sigphi*sigphi)?
			//get n subclusters
			int nSCs = jet.GetNConstituents();
			//if 1 subcluster, take dR to be 1 sigma
			if(nSCs < 2){
				Matrix cov, mu;
				jet.GetClusterParams(mu, cov);
				return max(sqrt(cov.at(0,0)), sqrt(cov.at(1,1)));

			}


			Jet subjet_i, subjet_j;			

			for(int i = 0; i < nSCs; i++){
				for(int j = i; j < nSCs; j++){
					//get center in eta, phi for subcluster i and j
					subjet_i = jet.GetConstituent(i);
					subjet_j = jet.GetConstituent(j);
					dr = subjet_i.deltaR(subjet_j);
					if(dr > maxDr) maxDr = dr;
				}
			}	
			return maxDr;
		}
	

		void CalcMMAvgPhiTime(BasePDFMixture* model, double& phi, double& t){
			int kmax = model->GetNClusters();
			phi = 0;
			t = 0;
			double pi, ws;
			map<string, Matrix> params;
			for(int k = 0; k < kmax; k++){
				params = model->GetPriorParameters(k);
				phi += params["pi"].at(0,0)*params["mean"].at(1,0);
				t += params["pi"].at(0,0)*params["mean"].at(2,0);
				ws += params["pi"].at(0,0);
			}
			phi /= ws;
			t /= ws; 
		}
		void MakeProcCats(string sample, bool leadsep = true){
			//total
			procCat tot(_hists1D, _hists2D);
			tot.ids = {-999};
			_procCats.push_back(tot);	
			
			if(sample.find("GMSB") != string::npos){
				//notSunm
				procCat notSunm(_hists1D, _hists2D, "notSunm","notSunm", leadsep);
				//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
				notSunm.ids = {97, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
				_procCats.push_back(notSunm);
				
				//signal
				procCat sig(_hists1D, _hists2D, "chiGam","#Chi^{0} #rightarrow #gamma", leadsep);
				sig.ids = {22};
				_procCats.push_back(sig);
			}
			else if(sample.find("JetHT") != string::npos){
				//data
				procCat jetht(_hists1D, _hists2D, "JetHT", "JetHT", leadsep);
				jetht.ids = {-999};
				_procCats.push_back(jetht);
			}
			else if(sample.find("GJets") != string::npos){
				//data
				procCat gjets(_hists1D, _hists2D, "GJets", "GJets", leadsep);
				gjets.ids = {-999};
				_procCats.push_back(gjets);
			}
			else if(sample.find("ttbar") != string::npos){
				procCat ttbar(_hists1D, _hists2D, "ttbar", "t#bar{t}",leadsep);
				ttbar.ids = {-999};
				_procCats.push_back(ttbar);
			}
			else if(sample.find("QCD") != string::npos){
				procCat qcd(_hists1D, _hists2D, "QCD", "QCD multi-jets",leadsep);
				qcd.ids = {-999};
				_procCats.push_back(qcd);
			}
			else return;

		}
		void WriteOutput(TFile* ofile){
			WriteStackHists(ofile);
			WriteHists(ofile);
			string name;
			ofile->cd();
			for(int i = 0; i < (int)graphs.size(); i++){
				//name = graphs[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDRGraph(graphs[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
				graphs[i]->Write();
			}
			ofile->Close();
		}	
		void WriteHists(TFile* ofile){
			string name;

			ofile->cd();
			//for condor skims, histograms need to be written because these
			//can be hadded together (TCanvases can't)
			for(int i = 0; i < (int)_hists1D.size(); i++){
				//name = hists1D[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDRHist(hists1D[i], cv, name, name, "a.u.");
				//write cv to file			
			//	cv->SaveAs((fname+"/"+name+".pdf").c_str());
				//cv->Write();
				if(_hists1D[i]->GetEntries() == 0) continue;
				_hists1D[i]->Write();
			}
			for(int i = 0; i < (int)_hists2D.size(); i++){
				//name = hists2D[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDR2DHist(hists2D[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
				if(_hists2D[i]->GetEntries() == 0) continue;
				_hists2D[i]->Write();
			}


		}
		
		void WriteStackHists(TFile* ofile){
			ofile->cd();
			string name, dirname, histname;
			//write 1D hists
			//variables
			int nhists = _procCats[0].hists1D[0].size();
			for(int i = 0; i < nhists; i++){
				name = _procCats[0].hists1D[0][i]->GetName();
				_procCats[0].hists1D[0][i]->Write();
				//dirname = name.substr(0,name.rfind("_"+trs[0].methodName));
			//cout << "i: " << i << " name " << name << " making dir " << name+"_stack" << endl;
				//TDirectory* dir = ofile->mkdir((name+"_stack").c_str());
				//dir->cd();
				//make process breakdown directory
				TDirectory *dir2 = ofile->mkdir((name+"_procStack").c_str());
				//cout << "  making dir " << dir2->GetName() << endl;
				dir2->cd();
				for(int p = 1; p < _procCats.size(); p++){
				//loop over processes
					histname = _procCats[p].hists1D[0][i]->GetName();
					if(_procCats[p].hists1D[0][i] == nullptr) continue;
					//cout << "    proc " << _procCats[p].plotName << " hist " << _procCats[p].hists1D[0][i]->GetName() << " " << _procCats[p].hists1D[0][i]->GetTitle() << " entries " << _procCats[p].hists1D[0][i]->GetEntries() << endl;			
					if(_procCats[p].hists1D[0][i]->GetEntries() == 0 && ((histname.find("sigma") == string::npos && histname.find("mean") == string::npos) && histname.find("profile") == string::npos)){ continue; }
					//cout << "  n hists " << trs[j].procCats[0].hists1D[0].size() << endl;
					//cout << "writing " << _procCats[p].hists1D[0][i]->GetName() << " " << _procCats[p].hists1D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;;
					_procCats[p].hists1D[0][i]->Write();

				}
				ofile->cd(); 
			}
			//cout << "2D hists" << endl;
			//write 2D hists
			nhists = _procCats[0].hists2D[0].size();
			for(int i = 0; i < nhists; i++){
				name = _procCats[0].hists2D[0][i]->GetName();
				histname = _procCats[0].hists2D[0][i]->GetName();
				//write total method histogram outside process directory
				if(_procCats[0].hists2D[0][i] == nullptr) continue;
				if(_procCats[0].hists2D[0][i]->GetEntries() == 0 && histname.find("sigma") == string::npos){ continue; }
				//check if data can be run
				//if(histname.find("recoGen") != string::npos && _data) continue;
				cout << "writing hist " << _procCats[0].hists2D[0][i]->GetName() << endl;
				_procCats[0].hists2D[0][i]->Write();
				//write method as directory within directory
				TDirectory *dir2 = ofile->mkdir((name+"_procStack").c_str());
				//cout << "  making dir " << dir2->GetName() << endl;
				dir2->cd();
				//cout << "writing " << _procCats[0].hists2D[0][i]->GetName() << " " << name << " " <<  _procCats[0].hists2D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;
				for(int p = 1; p < _procCats.size(); p++){
					//loop over processes
					if(_procCats[p].hists2D[0][i] == nullptr) continue;
					if(_procCats[p].hists2D[0][i]->GetEntries() == 0 && dirname.find("meanRecoGenDeltaT") == string::npos){ continue; }//cout << "Histogram for proc " << _plotName << " not filled." << endl; continue; }
					//check if data can be run
					histname = _procCats[p].hists2D[0][i]->GetName();
					//if(histname.find("recoGen") != string::npos && _data) continue;
					//cout << "  n hists " << _procCats[0].hists1D[0].size() << endl;
					_procCats[p].hists2D[0][i]->Write();
				} 
				ofile->cd(); 
			}
		}


		
		//comp time distribution
		TH1D* comptime = new TH1D("comptime","comptime",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime = new TGraph();
		
		//predicted jet plots
		//0
		TH1D* nClusters = new TH1D("nPredJets","nPredJets",10,0,10);
		//1
		TH1D* nSubClusters = new TH1D("predJet_nSubClusters","predJet_nSubClusters",10,0,10);
		//2
		TH1D* predJet_subClusterEnergy = new TH1D("predJet_subClusterEnergy","predJet_subClusterEnergy",20,0,500);
		//3
		TH1D* predJet_subClusterTimeCenter = new TH1D("predJet_subClusterTimeCenter","predJet_subClusterTimeCenter",25,0,15);
		//4
		TH1D* predJet_subClusterEtaCenter = new TH1D("predJet_subClusterEtaCenter","predJet_subClusterEtaCenter",25,-1.8,1.8);
		//5
		TH1D* predJet_subClusterPhiCenter = new TH1D("predJet_subClusterPhiCenter","predJet_subClusterPhiCenter",25,-0.1,6.3);
		//6
		TH1D* predJet_dR = new TH1D("predJet_dR","predJet_dR",50,0,2);
		//7
		TH1D* predJet_energy = new TH1D("predJet_energy","predJet_energy",50,0,400);
		//8
		TH1D* predJet_pt = new TH1D("predJet_pt","predJet_pt",50,0,300);
		//9
		TH1D* predJet_mass = new TH1D("predJet_mass","predJet_mass",50,0,200);
		//10 - resolution of difference of pt between reco and gen jets as a function of gen jet energy
		TH1D* jetGenE_sigmaDeltaPt_predGen = new TH1D("jetGenE_sigmaDeltaPt_predGen","jetGenE_sigmaDeltaPt_predGen",5,0,100);
		//11 - # pred jets - # gen jets
		TH1D* predGen_nJets = new TH1D("predGen_diffNJets","predGen_diffNJets",20,-10,10);
		//for subclusters
		//12 - eta sigma
		TH1D* predJet_subClusterEtaSigma = new TH1D("predJet_subClusterEtaSigma","predJet_subClusterEtaSigma",25,0.01,0.09);
		//13 - phi sigma
		TH1D* predJet_subClusterPhiSigma = new TH1D("predJet_subClusterPhiSigma","predJet_subClusterPhiSigma",25,0.01, 0.09);
		//14 - time sigma
		TH1D* predJet_subClusterTimeSigma = new TH1D("predJet_subClusterTimeSigma","predJet_subClusterTimeSigma",25,0.01,0.09);
		//15 - eta-phi covariance
		TH1D* predJet_subClusteretaPhiCov = new TH1D("predJet_subClusteretaPhiCov","predJet_subClusteretaPhiCov",25,-1,1);
		//16 - time-eta covariance
		TH1D* predJet_subClustertimeEtaCov = new TH1D("predJet_subClustertimeEtaCov","predJet_subClustertimeEtaCov",25,-1,1);
		//17 - time-phi covariance
		TH1D* predJet_subClustertimePhiCov = new TH1D("predJet_subClustertimePhiCov","predJet_subClustertimePhiCov",25,-1,1);
		//2D plots
		//0 - 2D histogram for recoGen pT resolution as a function of gen jet energy 
		TH2D* jetGenE_diffDeltaPt_predGen = new TH2D("jetGenE_diffDeltaPt_predGen","jetGenE_diffDeltaPt_predGen;jet_{gen} E (GeV);#Delta p_{T}_{pred, gen} (GeV)",5,0,100,50,-50,50);


		//reco jet plots
		//18
		TH1D* nRecoJets = new TH1D("nRecoJets","nRecoJets",10,0,10);
		//19
		TH1D* recoJet_dR = new TH1D("recoJet_dR","recoJet_dR",50,0,2);
		//20
		TH1D* recoJet_energy = new TH1D("recoJet_energy","recoJet_energy",50,0,150);
		//21
		TH1D* recoJet_pt = new TH1D("recoJet_pt","recoJet_pt",50,0,150);
		//22
		TH1D* recoJet_mass = new TH1D("recoJet_mass","recoJet_mass",50,0,150);
		//23 - resolution of difference of pt between reco and gen jets as a function of gen jet energy
		TH1D* jetGenE_sigmaDeltaPt_recoGen = new TH1D("jetGenE_sigmaDeltaPt_recoGen","jetGenE_sigmaDeltaPt_recoGen",5,0,100);
		//24 - # reco jets - # gen jets
		TH1D* recoGen_nJets = new TH1D("recoGen_diffNJets","recoGen_diffNJets",20,-10,10);
		//25 - reco jet pt/gen jet pt
		TH1D* recoGen_jetPtRatio = new TH1D("recoGen_jetPtRatio","recoGen_jetPtRatio",20,0,4);
		

		//2D plots
		//1 - 2D histogram for recoGen pT resolution as a function of gen jet energy 
		TH2D* jetGenE_diffDeltaPt_recoGen = new TH2D("jetGenE_diffDeltaPt_recoGen","jetGenE_diffDeltaPt_recoGen;jet_{gen} E (GeV);#Delta p_{T}_{reco, gen} (GeV)",5,0,100,50,-50,50);
		//2 - 2D histogram of gen pT vs reco pT
		TH2D* genPt_recoPt = new TH2D("genPt_recoPt","genPt_recoPt;genpt;recopt",50,5,50,50,5,50);

	
		void SetSmear(bool t){ _smear = t; }
		void SetTimeSmear(bool t){ _timesmear = t; }
		void SetOutfile(string fname){ _oname = fname; }
		void SetTransferFactor(double gev){
			_gev = gev;
			_prod->SetTransferFactor(_gev);
		}
		void SetEventRange(int evti, int evtj){ _evti = evti; _evtj = evtj; }
		void SetVerbosity(int verb){_verb = verb;}
		int _verb;

	private:
		string _oname;
		vector<TH1D*> _hists1D;
		vector<TH2D*> _hists2D;
		vector<TGraph*> graphs;
		vector<Jet> _phos; //photons for event
		vector<procCat> _procCats;
		vector<node*> _trees;
		vector<Jet> _predJets;
		vector<Jet> _genjets;
		vector<Jet> _recojets;
		bool _smear, _timesmear;
		enum Strategy{
			//Delauney strategy - NlnN time - for 2pi cylinder
			NlnN = 0,
			//traditional strategy - N^2 time
			N2 = 1,
		};
		//clustering strategy - N^2 or NlnN
		Strategy _strategy;
		ReducedBaseSim* _base = nullptr;
		int _nEvts;
		JetSimProducer* _prod = nullptr;
		bool _data;
		bool _debug;
		int _evti, _evtj;
		double _gev;
		double _c = 29.9792458; // speed of light in cm/ns
		double _radius; //radius of detector set by rhs in event (used for constructing jets)
		double _pvx, _pvy, _pvz;	


		double dR(double eta1, double phi1, double eta2, double phi2){
			return sqrt((eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2));
		}
		
};
#endif
