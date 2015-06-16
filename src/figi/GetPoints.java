package figi;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class GetPoints {

  public float[][] getTxf(int nt, float[][][] g) {
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
    float[] t = new float[nt];
    float[] x = new float[nt];
    float[] f = new float[nt];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int np = 0;
      int j1 = 0;
      for (int i1=0; i1<n1; ++i1) {
        if(g[i3][i2][i1]>0.0f) {np++;}
        if(g[i3][i2][i1]==0.0f) {j1++;}
      }
      if(np>nt) {
        int k = 0;
        for (; j1<nt+j1; j1++) {
          t[k] = k;
          x[k] = 172f;
          f[k] = g[i3][i2][j1];
          k++;
        }
      }
    }}
    return new float[][]{t,x,f};
  }

  public float[][] getCoordinates(float[][][] g) {
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
    int np = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(g[i3][i2][i1]>0.0f)
        np++;
    }}}
    int ip = 0;
    float[][] ks = new float[3][np];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(g[i3][i2][i1]>0.0f) {
        ks[0][ip] = i1;
        ks[1][ip] = i2;
        ks[2][ip] = i3;
        ip++;
      }
    }}}
    return ks;
  }

  public float[] getValues(float[][][] g) {
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
    int np = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(g[i3][i2][i1]>0.0f)
        np++;
    }}}
    int ip = 0;
    float[] fk = new float[np];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(g[i3][i2][i1]>0.0f) {
        fk[ip] = g[i3][i2][i1];
        ip++;
      }
    }}}
    return fk;
  }

  public float[][][] getWeights(float[][][] g) {
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
    float fnull = 0.0f;
    float[][][] ds = new float[n3][n2][n1];
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(fnull,g,ds,k1,k2,k3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(ds[i3][i2][i1]<=1.0f) {
        ds[i3][i2][i1] = 0.1f;
      }
    }}}
    System.out.println("pow="+log(5f)/log(max(ds)));
    return div(1,pow(ds,log(5f)/log(max(ds))));
    //return fillfloat(1f,n1,n2,n3);
  }

  public float[][] getWeights(Sampling s1, Sampling s2, float[] x1, float[] x2) 
  {
    int n2 = s2.getCount();
    int n1 = s1.getCount();
    float fnull = 0.0f;
    float[][] ds = new float[n2][n1];
    short[][] k1 = new short[n2][n1];
    short[][] k2 = new short[n2][n1];
    float[][] g  = new float[n2][n1];
    int np = x1.length;
    for (int ip=0; ip<np; ++ip) {
      int i1 = s1.indexOfNearest(x1[ip]);
      int i2 = s2.indexOfNearest(x2[ip]);
      g[i2][i1] = 1.0f;
    }
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(fnull,g,ds,k1,k2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(ds[i2][i1]<=1.0f) {
        ds[i2][i1] = 0.1f;
      } 
    }}
    System.out.println("pow="+log(6.67f)/log(max(ds)));
    return div(1,pow(ds,log(10f)/log(max(ds))));
    //return div(1,pow(ds,0.5f));
  }

}
