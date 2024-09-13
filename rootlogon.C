{
//gROOT->ProcessLine(" #include "\
//gStyle->SetPalette(kCandy);
gSystem->Load("lib/libvecDict.so");
gROOT->ProcessLine(".x SetPalette.C(2)");
gStyle->SetPadGridX(true);
gStyle->SetPadGridY(true);
}
