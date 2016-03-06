package mef;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import static mef.FaultGeometry.*;

/**
 * remove seismic reflections with slopes
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.03.06
 */

public class ReflectionRemove {

  public float[][][] apply(float[][][] p2, float[][][] p3, float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][][] g = new float[n3][n2][n1];
    float[] fm = new float[n1];
    float[] tm = new float[n1];
    float[] fp = new float[n1];
    float[] tp = new float[n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int i2m = (_bias== 0.5f || i2==0   )?i2:i2-1;
      int i2p = (_bias==-0.5f || i2==n2-1)?i2:i2+1;
      int i3m = (_bias== 0.5f || i3==0   )?i3:i3-1;
      int i3p = (_bias==-0.5f || i3==n3-1)?i3:i3+1;
      if (i2m<i2p && i3m<i3p) {
        float dx2 = 0.5f*(i2p-i2m);
        float dx3 = 0.5f*(i3p-i3m);
        for (int i1=0; i1<n1; ++i1) {
          tm[i1] = i1-p2[i3][i2][i1]*dx2-p3[i3][i2][i1]*dx3;
          tp[i1] = i1+p2[i3][i2][i1]*dx2+p3[i3][i2][i1]*dx3;
        }
        _si.interpolate(n1,1.0,0.0,f[i3m][i2m],n1,tm,fm);
        _si.interpolate(n1,1.0,0.0,f[i3p][i2p],n1,tp,fp);
        float scale = 0.5f/dx2;
        for (int i1=0; i1<n1; ++i1)
          g[i3][i2][i1] = scale*(fp[i1]-fm[i1]);
      }
    }}
    return g;
  }

  public float[][][] applyX(float[][][] p2, float[][][] p3, float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][][] g = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      applyFilter(p3[i3],f[i3],g[i3]);
    }
    for (int i2=0; i2<n2; ++i2) {
      float[][] f2i = new float[n3][n1];
      float[][] g2i = new float[n3][n1];
      float[][] p2i = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) {
        f2i[i3] = g[i3][i2];
        p2i[i3] = p2[i3][i2];
      }
      applyFilter(p2i,f2i,g2i);
      for (int i3=0; i3<n3; ++i3) {
        g[i3][i2] = g2i[i3];
      }
    }
    return g;
  }


  public void applyFilter(float[][] p, float[][] f, float[][] g) {
    int n1 = p[0].length;
    int n2 = p.length;
    float[] fm = new float[n1];
    float[] tm = new float[n1];
    float[] fp = new float[n1];
    float[] tp = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      int i2m = (_bias== 0.5f || i2==0   )?i2:i2-1;
      int i2p = (_bias==-0.5f || i2==n2-1)?i2:i2+1;
      if (i2m<i2p) {
        float dx2 = 0.5f*(i2p-i2m);
        for (int i1=0; i1<n1; ++i1) {
          tm[i1] = i1-p[i2][i1]*dx2;
          tp[i1] = i1+p[i2][i1]*dx2;
        }
        _si.interpolate(n1,1.0,0.0,f[i2m],n1,tm,fm);
        _si.interpolate(n1,1.0,0.0,f[i2p],n1,tp,fp);
        float scale = 0.5f/dx2;
        for (int i1=0; i1<n1; ++i1)
          g[i2][i1] = scale*(fp[i1]-fm[i1]);
      }
    }
  }

  private float _bias = 0.0f;
  private SincInterpolator _si = new SincInterpolator();


}

