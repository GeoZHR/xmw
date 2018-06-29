package tjxd;

import util.*;
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

  public float[][][] fill(float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] ff = copy(fx);
    int[][] bx1 = new int[n3][n2];
    int[][] tx1 = new int[n3][n2];
    float[][] fxb = new float[n3][n2];
    float[][] fxt = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int b1 = n1-1;
      for (int i1=n1-1; i1>=0; --i1) {
	float fxi = fx[i3][i2][i1];
        if(fxi!=0f) {
          bx1[i3][i2] = i1;
          fxb[i3][i2] = fxi;
	  break;
	}
      }
      int t1 = 0;
      for (int i1=0; i1<n1; ++i1) {
	float fxi = fx[i3][i2][i1];
        if(fxi!=0f) {
          tx1[i3][i2] = i1;
	  fxt[i3][i2] = fxi;
	  break;
	}
      }
    }}
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(4);
    rgf.apply00(fxt,fxt);
    rgf.apply00(fxb,fxb);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int b1 = bx1[i3][i2];
      int t1 = tx1[i3][i2];
      for (int i1=0; i1<=t1; i1++)
	ff[i3][i2][i1] = fxt[i3][i2];
      for (int i1=b1; i1<n1; i1++)
	ff[i3][i2][i1] = fxb[i3][i2];
    }}
    return ff;
  }

  public float[][][] fillLowVelocity(
    float i1min, float i1max, float vmin, float vmax, 
    float[][][] ex, float[][][] vx) {
    int n3 = vx.length;
    int n2 = vx[0].length;
    int n1 = vx[0][0].length;
    float[][][] vf = copy(vx);
    float dv = vmax-vmin;
    float d1 = i1max-i1min;
    dv /= d1;
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      if(ex[i3][i2][i1]>0f){
        vf[i3][i2][i1] = vmin+(i1-i1min)*dv;
      }
    }}}
    return vf;
  }

  public float[][][] fillLowVelocity(
    Sampling s1i, float i1min, float i1max, float vmin, float vmax, 
    float[][][] ex, float[][][] vx) {
    int n3 = vx.length;
    int n2 = vx[0].length;
    int n1 = vx[0][0].length;
    float[][][] vf = copy(vx);
    float dv = vmax-vmin;
    float d1 = i1max-i1min;
    dv /= d1;
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      double x1i = s1i.getValue(i1);
      int k1 = _s1.indexOfNearest(x1i);
      if(ex[i3][i2][k1]>0f){
        vf[i3][i2][i1] = vmin+(i1-i1min)*dv;
      }
    }}}
    return vf;
  }


  public void karstThreshold(float[][][] ex) {
    int n3 = ex.length;
    int n2 = ex[0].length;
    int n1 = ex[0][0].length;
    //zero the top (<450 samples)
    for (int i3=3; i3<n3-3; i3++) {
    for (int i2=3; i2<n2-3; i2++) {
    for (int i1=0; i1<490;  i1++) {
	    ex[i3][i2][i1] = 0.0f;
    }}}
    //zero the valued near boundaries
    for (int i3=0; i3<3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      ex[i3][i2][i1] = 0f;
    }}}
    for (int i3=n3-3; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      ex[i3][i2][i1] = 0f;
    }}}
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<3; i2++) {
    for (int i1=0; i1<n1; i1++) {
      ex[i3][i2][i1] = 0f;
    }}}
    for (int i3=0; i3<n3; i3++) {
    for (int i2=n2-3; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      ex[i3][i2][i1] = 0f;
    }}}
    for (int i3=3; i3<n3-3; i3++) {
    for (int i2=3; i2<n2-3; i2++) {
    for (int i1=490; i1<900;  i1++) {
      if(ex[i3][i2][i1]<0.18f) {
	      ex[i3][i2][i1] = 0.0f;
      }
    }}}

    for (int i3=3; i3<n3-3; i3++) {
    for (int i2=3; i2<n2-3; i2++) {
    for (int i1=900; i1<980;  i1++) {
      if(ex[i3][i2][i1]<0.21f) {
	      ex[i3][i2][i1] = 0.0f;
      }
    }}}

    for (int i3=3; i3<n3-3; i3++) {
    for (int i2=3; i2<n2-3; i2++) {
    for (int i1=980; i1<1160;  i1++) {
      if(ex[i3][i2][i1]<0.16f) {
	      ex[i3][i2][i1] = 0.0f;
      }
    }}}

    for (int i3=3; i3<n3-3; i3++) {
    for (int i2=3; i2<n2-3; i2++) {
    for (int i1=1160; i1<1280;  i1++) {
      if(ex[i3][i2][i1]<0.23f) {
	      ex[i3][i2][i1] = 0.0f;
      }
    }}}

    for (int i3=3; i3<n3-3; i3++) {
    for (int i2=3; i2<n2-3; i2++) {
    for (int i1=1280; i1<n1;  i1++) {
      if(ex[i3][i2][i1]<0.18f) {
	      ex[i3][i2][i1] = 0.0f;
      }
    }}}

  }

  public float[][][] resampling(final Sampling s1, final Sampling s1i, final float[][][] fx) {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int m1 = s1i.getCount();
    final float[][][] fi = new float[n3][n2][m1];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
	si.interpolate(s1,fx[i3][i2],s1i,fi[i3][i2]);
      }
    }});
    return fi;
  }

  public float[][] regridLogSamples(Sampling s1, float[] fx, float[] x1, float[] x2) {
    int n1 = s1.getCount();
    float dx1 = x1[1]-x1[0];
    float d1 = (float)s1.getDelta();
    int dm = round(d1/(dx1*2f));
    int[] ids = fillint(-1,n1);
    int nl = fx.length;
    for (int i1=0; i1<n1; i1++) {
      float x1i = (float)s1.getValue(i1);
      float[] dxs = abs(sub(x1,x1i));
      int[] id = new int[1];
      float dxm = min(dxs,id);
      ids[i1] = id[0];
    }
    float[] fxr = new float[n1];
    float[] x1r = new float[n1];
    float[] x2r = new float[n1];
    int j1 = 0;
    for (int i1=0; i1<n1; i1++) {
      int il = ids[i1];
      if(il<0) continue;
      x1r[j1] = x1[il];
      x2r[j1] = x2[il];
      fxr[j1] = fx[il];
      j1++;
    }
    fxr = copy(j1,0,fxr);
    x1r = copy(j1,0,x1r);
    x2r = copy(j1,0,x2r);
    return new float[][]{fxr,x1r,x2r};
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
  public float[][] apply(
    float[][] fx, float[][] x1, float[][] x2,  float[][] t) {
    final int n2 = t.length;
    final int n1 = t[0].length;
    final float[][] fi = new float[n2][n1];
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
          float fti = si.interpolate(_s1,_s2,t,x1i,x2i);
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
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x2i = (float)_s2.getValue(i2);
        float txi = t[i2][i1];
        float fxs = 0f;
        float scs = 0f;
        for (int iw=0; iw<nw; iw++) {
          int np = ft[iw].length;
          CubicInterpolator cii = cis.get(iw);
          if(txi>ft[iw][0]&&txi<ft[iw][np-1]) {
            float dx2 = x2i-x2[iw][0];
            float rxi = sqrt(dx2*dx2);
            float fxi = cii.interpolate(txi);
            float sci = radialBasis(rxi);
            scs += sci;
            fxs += fxi*sci;
          }
        }
        if(scs>0f) fi[i2][i1] = fxs/scs;
      }
    }});
    return fi;
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
    List<CubicInterpolator> ciw = new ArrayList<CubicInterpolator>();
    for (int iw=0; iw<nw; iw++) {
      int nl = x1[iw].length;
      float[] ws = new float[n1];
      for (int il=0; il<nl; ++il) {
        int i1 = _s1.indexOfNearest(x1[iw][il]);
        ws[i1] = 1f;
      }
      int i2 = _s2.indexOfNearest(x2[iw][0]);
      int i3 = _s3.indexOfNearest(x3[iw][0]);
      CubicInterpolator cii = new CubicInterpolator(t[i3][i2],ws);
      ciw.add(cii);
    }
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
          CubicInterpolator cwi = ciw.get(iw);
          if(txi>ft[iw][0]&&txi<ft[iw][np-1]) {
            float dx2 = x2i-x2[iw][0];
            float dx3 = x3i-x3[iw][0];
            float rxi = sqrt(dx2*dx2+dx3*dx3);
            float fxi = cii.interpolate(txi);
            float swi = cwi.interpolate(txi);
            float sci = swi*radialBasis(rxi);
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
   * Compute interpolation.
   * @param fx  2d array of known well-log property values
   * @param x1  2d array of vertical coordinates of the known well-log points
   * @param x2  2d array of inline coordinates of the known well-log points
   * @param x3  2d array of crossline coordinates of the known well-log points
   * @param t   3d array of RGT volume computed from a seismic image
   * @return    3d array of interpolated well-log property values.
   */
  public float[][][] apply(
    Sampling s1i, float[][] fx, float[][] x1, float[][] x2, 
    float[][] x3, float[][][] t) {
    final int n3 = t.length;
    final int n2 = t[0].length;
    final int n1 = t[0][0].length;
    final int m1 = s1i.getCount();
    final float[][][] fi = new float[n3][n2][m1];
    int nw = x1.length;
    // interpolate RGT values at the well-log points
    float[][] ft = new float[nw][]; //array of RGT values at well-log points
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    List<CubicInterpolator> cis = new ArrayList<CubicInterpolator>();
    List<CubicInterpolator> ciw = new ArrayList<CubicInterpolator>();
    for (int iw=0; iw<nw; iw++) {
      int nl = x1[iw].length;
      float[] ws = new float[n1];
      for (int il=0; il<nl; ++il) {
        int i1 = _s1.indexOfNearest(x1[iw][il]);
        ws[i1] = 1f;
      }
      int i2 = _s2.indexOfNearest(x2[iw][0]);
      int i3 = _s3.indexOfNearest(x3[iw][0]);
      CubicInterpolator cii = new CubicInterpolator(t[i3][i2],ws);
      ciw.add(cii);
    }
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
      for (int i1=0; i1<m1; ++i1) {
        float x2i = (float)_s2.getValue(i2);
        float x3i = (float)_s3.getValue(i3);
        float x1i = (float)s1i.getValue(i1);
        float txi = si.interpolate(n1,_s1.getDelta(),_s1.getFirst(),t[i3][i2],x1i);
        float fxs = 0f;
        float scs = 0f;
        for (int iw=0; iw<nw; iw++) {
          int np = ft[iw].length;
          CubicInterpolator cii = cis.get(iw);
          CubicInterpolator cwi = ciw.get(iw);
          if(txi>ft[iw][0]&&txi<ft[iw][np-1]) {
            float dx2 = x2i-x2[iw][0];
            float dx3 = x3i-x3[iw][0];
            float rxi = sqrt(dx2*dx2+dx3*dx3);
            float fxi = cii.interpolate(txi);
            float swi = cwi.interpolate(txi);
            float sci = swi*radialBasis(rxi);
            scs += sci;
            fxs += fxi*sci;
          }
        }
        if(scs>0f) fi[i3][i2][i1] = fxs/scs;
      }}
    }});
    return fi;
  }

  public float[][][] timeToDepthFunction(
    Sampling st, final float[][][] vt) {
    final int n3 = vt.length;
    final int n2 = vt[0].length;
    final int nt = vt[0][0].length;
    final float[][][] zt = new float[n3][n2][nt];
    final float dt = (float)st.getDelta();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        float[] vti = vt[i3][i2];
        float[] zti = zt[i3][i2];
        for (int i1=1; i1<nt; ++i1)
          zti[i1] = zti[i1-1]+1000f*vti[i1-1]*dt*0.5f;//meter
      }
    }});
    return zt;
  }

  public float[][][] timeToDepth(
    final float dz, final float[][][] zt, final float[][][] vt) {
    clean(0.001f,zt); // make sure zt is vertically mononic
    float zmax = max(zt);
    final int n3 = vt.length;
    final int n2 = vt[0].length;
    final int nz = round(zmax/dz);
    final float[][][] vz = new float[n3][n2][nz];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        float[] zti = zt[i3][i2];
        float[] vti = vt[i3][i2];
        float[] vzi = vz[i3][i2];
        CubicInterpolator ci = new CubicInterpolator(zti,vti);
        for (int i1=1; i1<nz; ++i1)
          vzi[i1] = ci.interpolate(i1*dz);
      }
    }});
    return vz;
  }

  public float[][][] fillTop(int h, int b1, float[][][] vx) {
    int n3 = vx.length;
    int n2 = vx[0].length;
    int n1 = vx[0][0].length;
    int d1 = h*round((float)(_s1.getFirst()/_s1.getDelta()));
    int m1 = n1+d1;
    float[][][] vf = new float[n3][n2][m1];
    float[][] vb = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      vb[i3][i2] = vx[i3][i2][b1];
    }}
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(5);
    rgf.apply00(vb,vb);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<d1+b1; ++i1) {
      vf[i3][i2][i1] = vb[i3][i2];
    }}}
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=b1; i1<n1; ++i1) {
      vf[i3][i2][d1+i1] = vx[i3][i2][i1];
    }}}
    return vf;
  }


  /**
   * Interpolation weigths defined by an inverse multiquadric function.
   * @param r distance.
   * @return value of radial basis function.
   */
  public float radialBasis(float r) {
    return 1f/(1f+_epsilon*_epsilon*r*r);
  }

  private void clean(float dmin, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      float fxi = fx[i3][i2][i1];
      float fxm = fx[i3][i2][i1-1];
      float fxt = fxm+dmin;
      if(fxi<fxt) fx[i3][i2][i1] = fxt;
    }}}

  }


  private float _epsilon = 1f;
  private Sampling _s1,_s2,_s3;
}
