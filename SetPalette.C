void SetPalette(){


	int k = 2;
	cout << "Using palette for " << k << " clusters." << endl;

	//number of gradients in the palette
	int nColors = 100;
	//set color palette
	vector<double> stops, red, green, blue;
	//number of end point colors
	int nMainColors = k;

	if(k < 2){
		cout << "Error: please give at least 2 clusters for palette creation." << endl;
		return;
	}

	//2 clusters
	if(k == 2){
		//[0,1] values are R, G or B/255.
		//first color = light blue
		red.push_back(0.52);
		green.push_back(0.79);
		blue.push_back(0.96);
		//second color = light pink
		red.push_back(0.89);
		green.push_back(0.52);
		blue.push_back(0.96);
		
		//where to switch colors
		stops.push_back(0.0);
		stops.push_back(0.5);
		stops.push_back(1.0);

	}

	//3 clusters
	if(k == 3){
		//[0,1] values are R, G or B/255.
		//first color = light blue
		red.push_back(0.52);
		green.push_back(0.79);
		blue.push_back(0.96);
		
		//second color = light pink
		red.push_back(0.89);
		green.push_back(0.52);
		blue.push_back(0.96);
		
		//third color = light green
		red.push_back(0.52);
		green.push_back(0.95);
		blue.push_back(0.79);
		
		//where to switch colors
		stops.push_back(0.0);
		stops.push_back(0.3333);
		stops.push_back(0.6666);
		stops.push_back(1.0);

	}


	Int_t palette[100];
	Int_t fi = TColor::CreateGradientColorTable(nMainColors,&stops[0],&red[0],&green[0],&blue[0],nColors);
	for (int i=0;i<nColors;i++) palette[i] = fi+i;
	
//	gStyle->SetPalette(100,palette);






}
