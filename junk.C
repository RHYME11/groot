Double_t *fg(double *x, double *par){


  // y = p0 + p1*x*x +p2/x

  
  return par[0] + par[1]*x[0]*x[0];
}
