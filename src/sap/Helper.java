package sap;

public class Helper {

  public float[][][] transpose(float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = new float[n3][n1][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gx[i3][i1][i2] = fx[i3][n2-i2-1][i1];
    }}}
    return gx;
  }

  public float[][][] transpose13(float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = new float[n1][n2][n3];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gx[i1][i2][i3] = fx[n3-i3-1][i2][i1];
    }}}
    return gx;
  }


}
