package ssi;

import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterpolator;
import static edu.mines.jtk.util.ArrayMath.*;

public final class FlattenHelper {

  public FlattenHelper(Sampling sz) {
    _sz = sz;
    _lz = sz.getLast();
    _nz = sz.getCount();
    _fz = sz.getFirst();
    _dz = sz.getDelta();
  }

  public void shiftsToRgt(float[][][] u) {
    int nx = u.length;
    int ny = u[0].length;
    int nz = u[0][0].length;
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        for (int iz=0; iz<nz; ++iz) {
          float ui = u[ix][iy][iz];
          float zi = (float)_fz+iz*(float)_dz;
          u[ix][iy][iz] = zi+ui*(float)_dz;
        }
      }
    }
  }

  public void rgtToShifts(float[][][] u) {
    int nx = u.length;
    int ny = u[0].length;
    int nz = u[0][0].length;
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        for (int iz=0; iz<nz; ++iz) {
          float ui = u[ix][iy][iz];
          float zi = (float)_fz+iz*(float)_dz;
          u[ix][iy][iz] = (ui-zi)/(float)_dz;
        }
      }
    }
  }

  public float[][][] rgtToHorizonVolume(final float[][][] u) {
    cleanRGT(u);
    final int nx = u.length;
    final int ny = u[0].length;
    final int nz = u[0][0].length;
    final float[][][] x = new float[nx][ny][nz];
    final InverseInterpolator ii = new InverseInterpolator(_sz,_sz);
    Parallel.loop(nx,new Parallel.LoopInt() {
    public void compute(int ix) {
      for (int iy=0; iy<ny; ++iy) 
        ii.invert(u[ix][iy],x[ix][iy]);
    }});
    return x;
  }

  public float[][][] flattenByRgt(final float[][][] u, final float[][][] f) {
    final int nx = u.length;
    final int ny = u[0].length;
    final int nz = u[0][0].length;
    final float[][][] g = new float[nx][ny][nz];
    final float[][][] x = rgtToHorizonVolume(u);
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(nx,new Parallel.LoopInt() {
    public void compute(int ix) {
      for (int iy=0; iy<ny; ++iy)
        si.interpolate(_nz,_dz,_fz,f[ix][iy],_nz,x[ix][iy],g[ix][iy]);
    }});
    /*
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        for (int iz=0; iz<nz; ++iz) {
          float xi = x[ix][iy][iz];
          //if(xi>=nz)  {g[ix][iy][iz]=Float.NaN;}
          //if(xi<0.0f) {g[ix][iy][iz]=Float.NaN;}
          if(xi>=nz)  {g[ix][iy][iz]=50.0f;}
          if(xi<0.0f) {g[ix][iy][iz]=50.0f;}
          if(iz>=1) {
            float xm = x[ix][iy][iz-1];
            //if(xi-xm<0.40f) {g[ix][iy][iz] = Float.NaN;}
            if(xi-xm<0.40f) {g[ix][iy][iz] = 50.0f;}
          }
        }
      }
    }
    */
    return g;
  }

  public void applyForHiatus(float[][] hv, float[][] fg) {
    int n2 = hv.length;
    int n1 = hv[0].length;
    float fgm = max(fg)*2.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = hv[i2][i1];
        if(xi>(double)_lz) {fg[i2][i1]=fgm;}
        if(xi<(double)_fz) {fg[i2][i1]=fgm;}
        if(i1>0) {
          float xm = hv[i2][i1-1];
          float xd = xi-xm;
          if(xd<(float)(0.3*_dz)) {fg[i2][i1] = fgm;}
        }
      }
    }
  }

  public float[][] horizonExtraction
    (float[][][] r, float[][] seed) 
  {
    int nx = r.length;
    int ny = r[0].length;
    int nz = r[0][0].length;
    float sz = seed[0][0];
    int sy = (int)seed[1][0];
    int sx = (int)seed[2][0];
    float[][] sf = fillfloat(sz,ny,nx);
    float[][][] h = rgtToHorizonVolume(r);
    final SincInterpolator si = new SincInterpolator();
    float sr = si.interpolate(nz,1.0,0.0,r[sx][sy],sz);
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        sf[ix][iy] = si.interpolate(nz,1.0,0.0,h[ix][iy],sr);
      }
    }
    return sf;
  }

  public int[][][] constraintsFromUnconformities(
    float[][][] p2, float[][][] p3, float[][][] ep, float[][][] unc) {
    if (unc==null) {return null;}
    int ns = unc.length;
    int n1 = p2[0][0].length;
    int[][][] uc = new int[3][ns][]; 
    for (int is=0; is<ns; ++is) {
      float[][] unci = unc[is];
      boolean fullSurf = checkFullUnc(n1,unci);
      float rd = 0.0f;
      float[][] rgti = null;
      if(fullSurf) {
        rgti = computeRgtOnUnconformity(p2,p3,ep,copy(unci));
        rd = min(rgti)+2.5f;
      }
      int np = checkValidPoints(n1,unci);
      uc[0][is] = new int[np];
      uc[1][is] = new int[np];
      uc[2][is] = new int[np];
      int nx = unci.length;
      int ny = unci[0].length;
      int ip = 0; 
      if(fullSurf) {
        for (int ix=0; ix<nx; ++ix) {
          for (int iy=0; iy<ny; ++iy) {
            float iz = unci[ix][iy];
            float ri = rgti[ix][iy];
            if (iz==0.0f||iz==n1-1.0f) {
              continue;
            } else if (ri>rd) {
              uc[2][is][ip] = ix;
              uc[1][is][ip] = iy;
              uc[0][is][ip] = round(iz);
              ip ++;
            }
          }
        }
      } else {
        for (int ix=0; ix<nx; ++ix) {
          for (int iy=0; iy<ny; ++iy) {
            float iz = unci[ix][iy];
            if(Float.isNaN(iz)){
              continue;
            } else {
              uc[2][is][ip] = ix;
              uc[1][is][ip] = iy;
              uc[0][is][ip] = round(iz);
              ip ++;
            }
          }
        }
      }
    }
    return uc;
  }

  public void setWeightsFromUnconformities(
    float[][][] wp, int[][][] uc, float fi, float fn) {
    int n1 = wp[0][0].length;
    int nc = uc[0].length;
    for (int ic=0; ic<nc; ++ic) {
      int np = uc[0][ic].length;
      for (int ip=0; ip<np; ++ip) {
        int i1 = uc[0][ic][ip];
        int i2 = uc[1][ic][ip];
        int i3 = uc[2][ic][ip];
        int i1m1 = i1-1;
        int i1p1 = i1+1;
        wp[i3][i2][i1] = fi; 
        if(i1m1>=0) {wp[i3][i2][i1m1] = fn;}
        if(i1p1<n1) {wp[i3][i2][i1p1] = fn;}
      }
    }
  }


  private int checkValidPoints(int n1, float[][] unc) {
    int np = 0;
    int nx = unc.length;
    int ny = unc[0].length;
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        float uci =  unc[ix][iy];
        if (uci!=0.0f) {np++;}
        if (uci!=n1-1.0f) {np++;}
        if (!Float.isNaN(uci)) {np++;}
      }
    }
    return np;
  }


  private boolean checkFullUnc(int n1, float[][] sf) {
    int n3 = sf.length;
    int n2 = sf[0].length;
    int sumN = 0;
    int sumP = n3*n2;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float sfi = sf[i3][i2];
        if(sfi == 0.0f) {sumN++;}
        if(sfi == n1-1.0f) {sumN++;}
        if(Float.isNaN(sfi)) {sumN++;}
      }
    }
    float perN = (float)sumN/(float)sumP;
    System.out.println("pern="+perN);
    if(perN<0.05f) {return true;}
    else           {return false;}
  }

  private float[][] computeRgtOnUnconformity(
    float[][][] p2, float[][][] p3, float[][][] ep, float[][] sf) {
    int nx = sf.length;
    int ny = sf[0].length;
    float[][] p2d = new float[nx][ny];
    float[][] p3d = new float[nx][ny];
    float[][] epa = new float[nx][ny];
    RgtOnUnconformity ru = new RgtOnUnconformity();
    ru.setSmoothings(8.0f,8.0f);
    ru.setCG(0.01f,200);
    checkSlopeAndWeight(p2,p3,ep,sf,p2d,p3d,epa);
    ru.rgtFromNormals(epa,p2d,p3d,sf,false);
    return sf;
  }

  private void checkSlopeAndWeight(
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][] sf,float[][] p2d, float[][] p3d, float[][] epd) 
  {
    int nx = p2.length;
    int ny = p2[0].length;
    int nz = p2[0][0].length;
    for (int ix=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        float sfi = sf[ix][iy];
        int iz = round(sfi);
        int zm2 = iz-2; if(zm2<0)    {zm2=0;}
        int zp2 = iz+2; if(zp2>nz-1) {zp2=nz-1;}
        int zm1 = iz-1; if(zm1<0)    {zm1=0;}
        int zp1 = iz+1; if(zp1>nz-1) {zp1=nz-1;}
        float p2m1 = p2[ix][iy][zm1];
        float p2p1 = p2[ix][iy][zp1];
        float p3m1 = p3[ix][iy][zm1];
        float p3p1 = p3[ix][iy][zp1];
        float epm1 = ep[ix][iy][zm1];
        float epp1 = ep[ix][iy][zp1];
        float p2m2 = p2[ix][iy][zm2];
        float p2p2 = p2[ix][iy][zp2];
        float p3m2 = p3[ix][iy][zm2];
        float p3p2 = p3[ix][iy][zp2];
        float epm2 = ep[ix][iy][zm2];
        float epp2 = ep[ix][iy][zp2];
        p2d[ix][iy] = 0.5f*(p2m2-p2p2+p2m1-p2p1);
        p3d[ix][iy] = 0.5f*(p3m2-p3p2+p3m1-p3p1);
        epd[ix][iy] = 0.25f*(epm1+epp1+epm2+epp2);
      }
    }
  }




  private void cleanRGT(float[][][] u) {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    float dMin = 0.001f*(float)_dz;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float ui = u[i3][i2][i1  ];
          float um = u[i3][i2][i1-1];
          float ut = um+dMin;
          if (ui<ut) {u[i3][i2][i1]=ut;}
        }
      }
    }
  }


  private static void cleanShifts(float[][][] r) {
    int n1 = r[0][0].length;
    int n2 = r[0].length;
    int n3 = r.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          if (r[i3][i2][i1]<=r[i3][i2][i1-1]-0.99f)
            r[i3][i2][i1] = r[i3][i2][i1-1]-0.99f;
        }
      }
    }
  }


  private final Sampling _sz; // vertical sampling
  private final int    _nz ; // number of samples
  private final double _fz ; // first samples
  private final double _lz ; // first samples
  private final double _dz ; // sample deltas
}
