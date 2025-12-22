void SetPalette(int k = 6){


	//int k = 2;
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

	//dark red
	red.push_back(132./255.);
	green.push_back(28./255.);
	blue.push_back(38./255.);
	
	//reddish pink
	red.push_back(186./255.);
	green.push_back(39./255.);
	blue.push_back(74./255.);
	
	//fifth color = pink
	red.push_back(227./255.);
	green.push_back(132./255.);
	blue.push_back(255./255.);
	//first color = light blue
	red.push_back(132./255.);
	green.push_back(201./255.);
	blue.push_back(244./255.);
	
	//teal
	red.push_back(178./255.);
	green.push_back(236./255.);
	blue.push_back(225./255.);

	//light light teal	
	red.push_back(238./255.);
	green.push_back(251./255.);
	blue.push_back(249./255.);

	/*
	//[0,1] values are R, G or B/255.
	//first color = light blue
	red.push_back(0.52);
	green.push_back(0.79);
	blue.push_back(0.96);
	//second color = light pink
	red.push_back(0.89);
	green.push_back(0.52);
	blue.push_back(0.96);
	
	//third color = light purple
	red.push_back(0.67);
	green.push_back(0.52);
	blue.push_back(0.95);

	//fifth color = light orange
	red.push_back(0.95);
	green.push_back(0.796);
	blue.push_back(0.52);

	//fourth color = light green
	red.push_back(0.52);
	green.push_back(0.95);
	blue.push_back(0.79);

	//sixth color = light yellow
	red.push_back(0.95);
	green.push_back(0.89);
	blue.push_back(0.52);
	*/	

	//where to switch colors
	stops.push_back(0.0);
	stops.push_back(0.16666);
	stops.push_back(0.33333);
	stops.push_back(0.5);
	stops.push_back(0.66666);
	stops.push_back(0.83333);
	stops.push_back(1.0);

	Int_t palette[100];
	Int_t fi = TColor::CreateGradientColorTable(nMainColors,&stops[0],&red[0],&green[0],&blue[0],nColors);
	for (int i=0;i<nColors;i++) palette[i] = fi+i;
	






}
