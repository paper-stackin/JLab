Double_t background(Double_t *x, Double_t *par) 
{

	Double_t bg_ampl = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];

    Double_t logX = log(x[0] + 0.2); // Добавляем 0.1 для избежания отрицательных значений
    Double_t logNorm = bg_ampl / ( (x[0] + 0.2) * sigma * sqrt(2 * M_PI)) * exp(-0.5 * pow((logX - mean) / sigma, 2));

    return logNorm;
}

// Define the peak function
Double_t peak(Double_t *x, Double_t *par) 
{
    // par[0] - amplitude of signal
    // par[1] - mean
    // par[2] - sigma
    // par[3] - radiation decay
    Double_t gauss = par[0] * TMath::Gaus(x[0], par[1], par[2]);
    Double_t tail = TMath::Exp(-par[3] * (x[0] - par[1]));

    if(x[0] >= par[1])  return gauss * tail;
    else                return gauss;
}

// Define the combined function
Double_t combined(Double_t *x, Double_t *par) 
{
    return background(x, par) + peak(x, &par[3]);
}