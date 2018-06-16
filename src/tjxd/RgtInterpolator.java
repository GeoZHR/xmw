package tjxd;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Well-log interpolation guided by relative geologic time.
 * Algorithm:
 * (1) Compute a relative geologic time (RGT) volume from a seismic image.
 * (2) Extrapolate well-log properties following constant RGT valumes.
 * (3) The interpolation is computed as the distance-weighted property 
 * values that are extrapolated from multiple well logs.
 * Mask is zero for all samples where labs is less than small*gabs.
 */
public class RgtInterpolator {

  /**
   * Constructs an RgtInterpolator.
   * @param s1 vertical sampling of the 3D RGT/seismic volume
   * @param s2 inline sampling of the 3D RGT/seismic volume
   * @param s3 crossline sampling of the 3D RGT/seismic volume
   * @param epsilon scale in the radial basis function
   */
  public RgtInterpolator(Sampling s1, Sampling s2, Sampling s3, float epsilon) {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _epsilon = epsilon;
  }


  /**
   * Compute interpolation.
   * @param fx  2d array of known well-log property values
   * @param x1  2d array of vertical coordinates of the known well-log points
   * @param x2  2d array of inline coordinates of the known well-log points
   * @param x3  2d array of crossline coordinates of the known well-log points
   * @param t   3d array of RGT volume computed from a seismic image
   * @return    3d array of interpolated well-log property values.
   */
  public float[][][] apply(
    float[][] fx, float[][] x1, float[][] x2, float[][] x3, float[][][] t) {
    final int n3 = t.length;
    final int n2 = t[0].length;
    final int n1 = t[0][0].length;
    final float[][][] fi = new float[n3][n2][n1];
    int nw = x1.length;
    // interpolate RGT values at the well-log points
    float[][] ft = new float[nw][]; //array of RGT values at well-log points
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    List<CubicInterpolator> cis = new ArrayList<CubicInterpolator>();
    for (int iw=0; iw<nw; iw++) {
      ArrayList<Float> fxs = new ArrayList<Float>();
      ArrayList<Float> fts = new ArrayList<Float>();
      int nl = x1[iw].length;
      ft[iw] = new float[nl];
      float ftp = -n1;
      for (int il=0; il<nl; ++il) {
        float x1i = x1[iw][il];
        if(x1i>_s1.getFirst()&&x1i<_s1.getLast()) {
          float x2i = x2[iw][il];
          float x3i = x3[iw][il];
          float fti = si.interpolate(_s1,_s2,_s3,t,x1i,x2i,x3i);
          if(fti>ftp) {
            fts.add(fti);
            fxs.add(fx[iw][il]);
            ftp = fti;
          }
        }
      }
      int np = fts.size();
      System.out.println("np="+np);
      float[] fta = new float[np];
      float[] fxa = new float[np];
      for (int ip=0; ip<np; ++ip) {
        fta[ip] = fts.get(ip);
        fxa[ip] = fxs.get(ip);
      }
      ft[iw] = fta;
      CubicInterpolator cii = new CubicInterpolator(fta,fxa);
      cis.add(cii);
    }
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x2i = (float)_s2.getValue(i2);
        float x3i = (float)_s3.getValue(i3);
        float txi = t[i3][i2][i1];
        float fxs = 0f;
        float scs = 0f;
        for (int iw=0; iw<nw; iw++) {
          int np = ft[iw].length;
          CubicInterpolator cii = cis.get(iw);
          if(txi>ft[iw][0]&&txi<ft[iw][np-1]) {
            float dx2 = x2i-x2[iw][0];
            float dx3 = x3i-x3[iw][0];
            float rxi = sqrt(dx2*dx2+dx3*dx3);
            float fxi = cii.interpolate(txi);
            float sci = radialBasis(rxi);
            scs += sci;
            fxs += fxi*sci;
          }
        }
        if(scs>0f) fi[i3][i2][i1] = fxs/scs;
      }}
    }});
    return fi;
  }

  /**
   * Interpolation weigths defined by an inverse multiquadric function.
   * @param r distance.
   * @return value of radial basis function.
   */
  public float radialBasis(float r) {
    return 1f/(1f+_epsilon*_epsilon*r*r);
  }


  private float _epsilon = 1f;
  private Sampling _s1,_s2,_s3;
}
