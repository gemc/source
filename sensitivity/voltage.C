// ROOT function to try values of parameters


double DGauss(double *x, double *par)
{
	double stepTime = par[6];
	double Edep     = par[5];

	double t0   = par[0] + stepTime;     // start time of signal so that peak is t0 + rise
	double rise = par[1]/3;    // rise time, equal to sigma of first gaussian.
	double fall = par[2]/3;    // fall time, equal to sigma of second gaussian
	double ampl = Edep*par[3]/2;    // amplitude of the signal mV/MeV
	double peds = par[4];      // pedestal
	
	
	double peak = t0 + 3*rise;
	
	return peds - ampl*exp(-0.5*pow((x[0]-peak)/rise, 2)) - ampl*exp(-0.5*pow((x[0]-peak)/fall, 2));
}


void v()
{

	TF1 *signal = new TF1("DGauss", DGauss, 0, 800, 7);
	signal->SetParameter(0, 50);
	signal->SetParameter(1, 1);
	signal->SetParameter(2, 2);
	signal->SetParameter(3, 100);
	signal->SetParameter(4, -20);
	
	signal->SetParameter(5, 22.374);
	signal->SetParameter(6, 5.55);
	signal->SetNpx(10000);
	
	
	signal->Draw();
	
}


