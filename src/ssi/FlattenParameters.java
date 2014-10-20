package ssi;

public enum FlattenParameters {
  INSTANCE();

  // Smoothing parameters
  public double sigma1;
  public double sigma2;
  public double sigma3;
  public double wPower;
  public double slpMax;


  // horizons and unconformities for control
  public String[] hvName;
  public String[] ucName;
  
  private FlattenParameters() {
    // Set default values, these public fields are updated from the UI.
    this.hvName = null;
    this.ucName = null;
    this.sigma1 = 2.0;
    this.sigma2 = 1.0;
    this.sigma3 = 1.0;
    this.slpMax = 4.0;
    this.wPower = 12.0;
  }
}
