{

	TF2 *cosmic = new TF2("cosmic","pow([0], [1]*cos(x*3.141592/180))/([2]*y*y)", 0, 90, 1, 10);

	cosmic->SetParameter(0, 55.6);
	cosmic->SetParameter(1, 1.04);
	cosmic->SetParameter(2, 64);

	cosmic->Draw("surf2");
	
}