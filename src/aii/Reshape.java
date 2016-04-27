package aii;

public class Reshape {

  public Reshape(int j2, int j3, int n2, int n3) {
    _j2 = j2;
    _j3 = j3;
    _n2 = n2;
    _n3 = n3;
  }

  public float[][][] apply(float[][] hd, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][][] gx = new float[_n3][_n2][n1];
    for (int k2=0; k2<n2; ++k2) {
      int i3 = (int)hd[k2][1]-_j3;
      int i2 = (int)hd[k2][2]-_j2;
      gx[i3][i2] = fx[k2];
    }
    return gx;
  }

  public float[][][] resample(int d1, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    int m1 = (n1-1)/8+1;
    float[][][] gx = new float[n3][n2][m1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k1=0;
    for (int i1=0; i1<n1; i1+=8) {
      gx[i3][i2][k1] = fx[i3][i2][i1];
      k1++;
    }}}
    return gx;
  }
  private int _j2, _j3, _n2, _n3;
}
