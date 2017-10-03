package gu;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Kaust edge detection
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.08.11
 */
public class KaustEdge {

  public float[][][] directionalDifference(
    float[][][] gx, float[][][] v1, float[][][] v2, float[][][] v3,
    float[][][] w1, float[][][] w2, float[][][] w3) {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final float[][][] dd = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    final float[][][] g1 = new float[n3][n2][n1];
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    rgf.apply100(gx,g1);
    rgf.apply010(gx,g2);
    rgf.apply001(gx,g3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float v1i = v1[i3][i2][i1];
      float v2i = v2[i3][i2][i1];
      float v3i = v3[i3][i2][i1];
      float w1i = w1[i3][i2][i1];
      float w2i = w2[i3][i2][i1];
      float w3i = w3[i3][i2][i1];
      float g1i = g1[i3][i2][i1];
      float g2i = g2[i3][i2][i1];
      float g3i = g3[i3][i2][i1];
      float gvi = v1i*g1i+v2i*g2i+v3i*g3i;
      float gwi = w1i*g1i+w2i*g2i+w3i*g3i;
      dd[i3][i2][i1] = sqrt(gvi*gvi+gwi*gwi);
    }}
    }});
    return dd;
  }
}
