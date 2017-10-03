package gu;

<<<<<<< HEAD
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;


/**
 * Horizon display
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.08.11
 */
public class HorizonDisplay {
 

  public float[][][] heightRgb(
    ColorMap mp, float[][] sf) {
    int n3 = sf.length;
    int n2 = sf[0].length;
    float[] sa = new float[n3*n2];
    int k = 0;
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      sa[k++] = sf[i3][i2];
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
  }

  public float[][][] amplitudeRgb(
    ColorMap mp, float[][][] fx, float[][] sf) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[] sa = new float[n3*n2];
    int k = 0;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      sa[k++] = si.interpolate(n1,1,0,fx[i3][i2],sf[i3][i2]);
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
  }

=======
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Kaust edge detection
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.08.11
 */
public class KaustEdge {

  float[][][] directionalDifference(
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
>>>>>>> 31efcfe16c653f77ffe4ba3b9436970550efe51b
}
