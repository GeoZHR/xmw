package swt;

import java.util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Seismic well ties.
 * 
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.06.09
 */

public class SeismicWellTie {


  public float[][][] logsToArray(WellLog[] logs, float[] zs) {
    int nc = 2;
    float dz = 0.0f;
    int nl = logs.length;
    float fk = 0.0003048f; // 1 ft = 0.0003048 km
    float[] fzs = new float[nl];
    float[] lzs = new float[nl];
    for (int il=0; il<nl; ++il){
      float[] z = logs[il].z;
      float[] d = logs[il].d;
      float[] v = logs[il].v;
      int nz = z.length;
      int jz = 0;
      int lz = nz-1;
      while (jz<nz && d[jz]==_vnull && v[jz]==_vnull)
        ++jz;
      while (lz>=0 && d[jz]==_vnull && v[jz]==_vnull)
        --lz;
      fzs[il] = z[jz];
      lzs[il] = z[lz];
      dz = z[1]-z[0];
    }
    float fz = min(fzs);
    float lz = max(lzs);
    int nz = round((lz-fz)/dz)+1;
    float[][][] wa = fillfloat(_vnull,nz,nl,nc);
    for (int il=0; il<nl; ++il) {
      float[] z = logs[il].z;
      float[] d = logs[il].d;
      float[] v = logs[il].v;
      int jz = 0;
      while (jz<z.length && d[jz]==_vnull && v[jz]==_vnull)
        ++jz;
      int jw = round((z[jz]-fz)/dz);
      for( ; jz<z.length; jz++) {
        wa[0][il][jw] = d[jz];
        wa[1][il][jw] = v[jz];
        jw ++;
      }
    }
    dz *= fk;
    fz *= fk;
    zs[0] = nz;
    zs[1] = dz;
    zs[2] = fz;
    return wa;
  }

  public float[][] computeSyns(
    boolean simple, WellLog[] logs, float[][] ndft, float[][] ndfz) 
  {
    float fpeak = 35f;
    float q = 100.0f;

    int il = 0;
    float ds = 0.002f;
    int nl = logs.length;
    float[][] sa = new float[nl][];
    for (WellLog log:logs) {
      SynSeis.Model md = SynSeis.getModel(log);
      float fsi = md.tmin();
      float lsi = md.tmax();
      fsi = round(fsi/ds)*ds;
      lsi = round(lsi/ds)*ds;
      int nsi = round((lsi-fsi)/ds)+1;
      Sampling ssi = new Sampling(nsi,ds,fsi);
      Sampling szi = md.sz;
      ndft[il][0] = nsi;
      ndft[il][1] = ds;
      ndft[il][2] = fsi;
      ndfz[il][0] = szi.getCount();
      ndfz[il][1] = (float)szi.getDelta();
      ndfz[il][2] = (float)szi.getFirst();
      if (simple) {
        sa[il] = SynSeis.makeSimpleSeismogram(md,fpeak,ssi);
      } else {
        sa[il] = SynSeis.makeBetterSeismogram(md,q,fpeak,ssi);
      }
      il++;
    }
    return sa;
  }

  public float[][] getTimeDepth(WellLog[] logs) {
    int il = 0;
    int nl = logs.length;
    float[][] td = new float[nl][];
    for (WellLog log:logs) {
      SynSeis.Model md = SynSeis.getModel(log);
      td[il] = md.t;
      il++;
    }
    return td;
  }

  public float[][] applySynsFlat(
    Sampling sz, Sampling sm[], Sampling[] st, float[][] ndfw, 
    float[][] ss, float[][] tz, float[][] gt) 
  {
    int nl = sm.length;
    float[][] gw = new float[nl][];
    for (int il=0; il<nl; ++il) {
      int nz    = sz.getCount();
      double fz = sz.getFirst();
      double dz = sz.getDelta();
      double lm = sm[il].getLast();
      double fm = sm[il].getFirst();
      double dm = sm[il].getDelta();
      int nw = 0;
      double fw = 0.0;
      ArrayList<Float> gl = new ArrayList<Float>();
      SincInterpolator sit = new SincInterpolator();
      SincInterpolator sig = new SincInterpolator();
      for (int ik=0; ik<nz; ++ik) {
        float s = ss[il][ik]; 
        double w = ik*dz+fz;
        double z = w-s*dz;
        if(s!=_vnull && z>fm+dm && z<lm-2*dm) {
          if(nw==0){fw=w;}
          float t = sit.interpolate(sm[il],tz[il],z);
          float g = sig.interpolate(st[il],gt[il],t);
          gl.add(g);
          nw++;
        }
      }
      gw[il] = new float[nw];
      for (int iw=0; iw<nw; ++iw) {
        gw[il][iw] = gl.get(iw);
      }
      ndfw[il][0] = nw;
      ndfw[il][1] = (float)dz;
      ndfw[il][2] = (float)fw;
    }
    return gw;
  }

  private float _vnull = -999.2500f;

}

