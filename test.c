Double_t Bateman(Double_t *x, Double_t *par){
  // par[0] = generation of one isotope in the decay chain (parent = 1, daughter = 2, granddaughter = 3 ...
  // par[1] = partilce[gen] activity
  // par[i+1] = half-life of particle[i] (gen >= i >= 1);

  int gen = par[0];
  double mul = 1;
  double sum = 0;
  for(int i=1;i<=gen;i++){
    double lami = 0.693/par[i+1];
    if(i>=2 && (i<=gen-1)) mul = mul*lami;
    double denominator = 1;
    for(int j=1;j<=gen;j++){
      if(j==i) continue;
      double lamj = 0.693/par[j+1];
      denominator = denominator * (lamj - lami);
    }
    sum += TMath::Exp(-lami*x[0])/denominator;
  }
  if(gen==1) return par[1]*mul*sum;
  return (0.693/par[gen+1])*(par[1]*mul*sum);
}

Double_t DecayChain(double *x, double *par){
  // par[0] = generation of isotopes in the decay chain (parent = 1, daughter = 2 ...);
  // par[1] = activity;
  // par[i+2] = half-life of particle_i(i>0)

  int gen = par[0];
  std::vector<double *> vec_par;
  for(int i=1;i<=gen;i++){
    double subpar[i+2];
    subpar[0] = (double) i;
    for(int j=1;j<=i+1;j++){
      subpar[j] = par[j];
    }
    vec_par.push_back(subpar);
  }

  double sum = 0;
  for(int i=1;i<=gen;i++){
    sum += Bateman(x,vec_par[i-1]);
  }

  return sum;
}
