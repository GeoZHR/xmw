package swt;

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
    float NULL_VALUE = -999.25f;
    float[] fzs = new float[nl];
    float[] lzs = new float[nl];
    for (int il=0; il<nl; ++il){
      float[] z = logs[il].z;
      float[] d = logs[il].d;
      float[] v = logs[il].v;
      int nz = z.length;
      int jz = 0;
      int lz = nz-1;
      while (jz<nz && d[jz]==NULL_VALUE && v[jz]==NULL_VALUE)
        ++jz;
      while (lz>=0 && d[jz]==NULL_VALUE && v[jz]==NULL_VALUE)
        --lz;
      fzs[il] = z[jz];
      lzs[il] = z[lz];
      dz = z[1]-z[0];
    }
    float fz = min(fzs);
    float lz = max(lzs);
    int nz = round((lz-fz)/dz)+1;
    float[][][] wa = fillfloat(NULL_VALUE,nz,nl,nc);
    for (int il=0; il<nl; ++il) {
      float[] z = logs[il].z;
      float[] d = logs[il].d;
      float[] v = logs[il].v;
      int jz = 0;
      while (jz<z.length && d[jz]==NULL_VALUE && v[jz]==NULL_VALUE)
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

}

