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
  private int _j2, _j3, _n2, _n3;
}
