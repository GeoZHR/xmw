package slt;

import java.util.*;


import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;
import ipfx.FaultCell;
import static ipfx.FaultGeometry.*;


/**
 * Extract salt boundary surfaces from oriented salt points. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.10
 */

public class SaltSurfer {


     // Uses fault images to find cells, oriented points located on ridges.
  public FaultCell[] findCells( float fmin,
    float[][][] f, float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    // Vertical image boundaries are discontinuities that may look like
    // faults. If a fault appears to be near and nearly parallel to image
    // boundaries, then assume it is a boundary artifact and not truly a
    // fault.
    int imax = 5; // max number of samples considered to be near boundary
    float wwmax = 0.75f; // cosine of 30 degrees, squared

    // Loop over all samples. Construct cells for samples nearest to ridges.
    ArrayList<FaultCell> cellList = new ArrayList<FaultCell>();
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmi = f[i3m][i2 ];
        float[] fim = f[i3 ][i2m];
        float[] fip = f[i3 ][i2p];
        float[] fpi = f[i3p][i2 ];
        float[] fmm = f[i3m][i2m];
        float[] fpp = f[i3p][i2p];
        float[] fmp = f[i3m][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fii = f[i3 ][i2 ];
        float[] u1i = u1[i3][i2];
        float[] u2i = u2[i3][i2];
        float[] u3i = u3[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          float fmii = fmi[i1 ];
          float fimi = fim[i1 ];
          float fipi = fip[i1 ];
          float fpii = fpi[i1 ];
          float fmmi = fmm[i1 ];
          float fppi = fpp[i1 ];
          float fmpi = fmp[i1 ];
          float fpmi = fpm[i1 ];
          float fiii = fii[i1 ];
          float u1ii = u1i[i1 ];
          float u2ii = u2i[i1 ];
          float u3ii = u3i[i1 ];
          if (u2ii==0f&&u3ii==0f){continue;}
          if (u1ii>0.0f) {
            u1ii = -u1ii;
            u2ii = -u2ii;
            u3ii = -u3ii;
          }
          float piii = faultStrikeFromNormalVector(u1ii,u2ii,u3ii);
          // Most image samples will not have a fault cell.
          FaultCell cell = null;

          // Accumulators for ridge likelihoods and locations. Depending on
          // the limits on fault strike used below, we may find more than one
          // ridge.
          float nr = 0;
          float fl = 0.0f;
          float d2 = 0.0f;
          float d3 = 0.0f;

          // If S-N ridge, ...
          if ((fipi<fiii && fimi<fiii) &&
              ((337.5f<=piii || piii<= 22.5f) || 
               (157.5f<=piii && piii<=202.5f))) {
            float f1 = 0.5f*(fipi-fimi); // 1st derivative
            float f2 = fipi-2.0f*fiii+fimi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=fmin) {
              if (imax<=i2 && i2<n2-imax || u2ii*u2ii<=wwmax) {
                fl += fr;
                d2 += dr;
                nr += 1;
              }
            }
          }

          // If SW-NE ridge, ...
          if ((fmpi<fiii && fpmi<fiii) &&
              (( 22.5f<=piii && piii<= 67.5f) || 
               (202.5f<=piii && piii<=247.5f))) {
            float f1 = 0.5f*(fmpi-fpmi); // 1st derivative
            float f2 = fmpi-2.0f*fiii+fpmi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=fmin) {
              if ((imax<=i2 && i2<n2-imax || u2ii*u2ii<=wwmax) &&
                  (imax<=i3 && i3<n3-imax || u3ii*u3ii<=wwmax)) {
                fl += fr;
                d2 += dr;
                d3 -= dr;
                nr += 1;
              }
            }
          }

          // If W-E ridge, ...
          if ((fpii<fiii && fmii<fiii) &&
              (( 67.5f<=piii && piii<=112.5f) ||
               (247.5f<=piii && piii<=292.5f))) {
            float f1 = 0.5f*(fpii-fmii); // 1st derivative
            float f2 = fmii-2.0f*fiii+fpii; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=fmin) {
              if (imax<=i3 && i3<n3-imax || u3ii*u3ii<=wwmax) {
                fl += fr;
                d3 += dr;
                nr += 1;
              }
            }
          }

          // If NW-SE ridge, ...
          if ((fppi<fiii && fmmi<fiii) &&
              ((112.5f<=piii && piii<=157.5f) || 
               (292.5f<=piii && piii<=337.5f))) {
            float f1 = 0.5f*(fppi-fmmi); // 1st derivative
            float f2 = fppi-2.0f*fiii+fmmi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=fmin) {
              if ((imax<=i2 && i2<n2-imax || u2ii*u2ii<=wwmax) &&
                  (imax<=i3 && i3<n3-imax || u3ii*u3ii<=wwmax)) {
                fl += fr;
                d2 += dr;
                d3 += dr;
                nr += 1;
              }
            }
          }

          // If at least one ridge, construct a cell and add to list.
          if (nr>0) {
            fl /= nr;
            d2 /= nr;
            d3 /= nr;
            float tiii = faultDipFromNormalVector(u1ii,u2ii,u3ii);
            cell = new FaultCell(i1,i2+d2,i3+d3,fl,piii,tiii);
            cellList.add(cell);
          }
        }
      }
    }
    return cellList.toArray(new FaultCell[0]);
  }

  public float[][][] scalarField(
    int n1, int n2, int n3, FaultCell[] cells) {
    float[][][] fx = new float[n3][n2][n1];
    for (FaultCell cell:cells) {
      if(cell==null) {continue;}
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      fx[i3][i2][i1] = cell.getFl();
    }
    return fx;
  }

  public FaultCell[] findPoints (
    final float fmin,
    final float[][][] fx, final float[][][] u1, 
    final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final float[][][] fm = new float[n3][n2][n1];
    final float[][][] fp = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    final ArrayList<FaultCell> cellList = new ArrayList<FaultCell>();
    Parallel.loop(1,n3-1,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=1; i2<n2-1 ;++i2) {
    for (int i1=1; i1<n1-1 ;++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float x1m = i1-u1i;
      float x2m = i2-u2i;
      float x3m = i3-u3i;
      float x1p = i1+u1i;
      float x2p = i2+u2i;
      float x3p = i3+u3i;
      fm[i3][i2][i1] = si.interpolate(s1,s2,s3,fx,x1m,x2m,x3m);
      fp[i3][i2][i1] = si.interpolate(s1,s2,s3,fx,x1p,x2p,x3p);
    }}}});
    for (int i3=1; i3<n3-1; ++i3) {
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float fxi = fx[i3][i2][i1];
      float fxm = fm[i3][i2][i1];
      float fxp = fp[i3][i2][i1];
      if (fxi>fxm && fxi>fxp && fxi>fmin && abs(u1i)<0.9f) {
        if (u1i>0.0f) {
          u1i = -u1i;
          u2i = -u2i;
          u3i = -u3i;
        }
        float t = faultDipFromNormalVector(u1i,u2i,u3i);
        float p = faultStrikeFromNormalVector(u1i,u2i,u3i);
        FaultCell cell = new FaultCell(i1,i2,i3,fxi,p,t);
        cellList.add(cell);
      }
    }}}
    return cellList.toArray(new FaultCell[0]);
  }

}
