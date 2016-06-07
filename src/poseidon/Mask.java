package poseidon;

public class Mask {
  
  public void setValue(float fv, float[][][] gx, float[][][] fx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float gxi = gx[i3][i2][i1];
      if (gxi==0f) {
        fx[i3][i2][i1] = fv;
      }
    }}}
  }

}
