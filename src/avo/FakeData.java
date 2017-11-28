/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package avo;
import java.util.Random;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import util.*;

import static ipf.FaultGeometry.*;

/**
 * Generates fake seismic data with faults and horizons.
 * <em>
 * Jacobians of functions used in folding and faulting have been implemented
 * but not tested. Therefore, beware of errors in calculated slopes p2 and p3.
 * </em>
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.05.05
 */
public class FakeData {

  /**
   * Plots specified fake data. Data type is specified by method name. For
   * example, type "seismicAndSlopes2d2014A" corresponds to the method
   * seismicAndSlopes2d2014A().
   */
  public static void main(final String[] args) {
    SwingUtilities.invokeLater(new Runnable(){
      public void run() {
        go(args);
      }
    });
  }
  private static void go(String[] args) {
    if (args.length==0 || args[0].equals("seismicAndSlopes2d2014A")) {
      float[][][] gp = seismicAndSlopes2d2014A(0.0);
      float[][] g = gp[0];
      float[][] p = gp[1];
      SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
      PixelsView pv = sp.addPixels(g);
      //pv.setInterpolation(PixelsView.Interpolation.LINEAR);
      sp.setSize(700,700);
    } else if (args[0].equals("seismicAndSlopes3d2014A")) {
      String[] sequences = {"OA"};
      int[] nplanars = {0};
      boolean[] conjugates = {false};
      boolean[] conicals = {true};
      boolean[] impedances = {false};
      boolean[] wavelets = {true};
      double[] noises = {0.5};
      int n = sequences.length;
      for (int i=0; i<n; ++i) {
        float[][][][] f = seismicAndSlopes3d2014A(
            sequences[i],nplanars[i],conjugates[i],conicals[i],
            impedances[i],wavelets[i],noises[i]);
        trace(" f min="+min(f[0])+" max="+max(f[0]));
        trace("p2 min="+min(f[1])+" max="+max(f[1]));
        trace("p3 min="+min(f[2])+" max="+max(f[2]));
        SimpleFrame frame = new SimpleFrame();
        frame.setSize(700,700);
        ImagePanelGroup ipg = frame.addImagePanels(f[0]);
        ipg.setSlices(100,51,24);
        if (impedances[i]) {
          ipg.setColorModel(ColorMap.JET);
        } else {
          ipg.setColorModel(ColorMap.GRAY);
        }
        //ViewCanvas vc = frame.getViewCanvas();
        OrbitView ov = frame.getOrbitView();
        //vc.setBackground(Color.WHITE);
        ov.setAzimuthAndElevation(40.0,25.0);
        ov.setScale(2.2);
        ov.setTranslate(new Vector3(-0.0102,-0.0508,0.0395));
      }
    } else {
      System.out.println("unrecognized type of fake data");
    }
  }

  /**
   * Returns a fake 3D seismic image with slopes.
   * @param noise rms of noise (relative to signal) added to the image.
   * @return array of arrays {f,p2,p3} with image f and slopes p2 and p3.
   */
  public static float[][][][] seismicAndSlopes3d2014A(double noise) {
    return seismicAndSlopes3d2014A("OA",3,false,false,false,true,noise);
  }

  /**
   * Returns a fake 3D seismic image with slopes.
   * @param sequence string of 'O' and 'A' for fOlding and fAulting.
   * @param nplanar number of planar faults.
   * @param conjugate true, to make one of the faults a conjugate fault.
   * @param conical true, to include a conical fault.
   * @param impedance true, for impedance instead of reflectivity.
   * @param wavelet true, for wavelet; false, for no wavelet.
   * @param noise rms of noise (relative to signal) added to the image.
   * @return array of arrays {f,p2,p3} with image f and slopes p2 and p3.
   */
  public static float[][][][] seismicAndSlopes3d2014A(
      String sequence, int nplanar, boolean conjugate, boolean conical,
      boolean impedance, boolean wavelet, double noise) {
    int n1 = 101;
    int n2 = 152;
    int n3 = 153;
    int m1 = n1+50;

    // Number of episodes in deformation sequence of folding and faulting.
    int ns = sequence.length();
    int no = 0;
    int na = 0;
    for (int js=0; js<ns; ++js) {
      if (sequence.charAt(js)=='O') {
        ++no;
      } else if (sequence.charAt(js)=='A') {
        ++na;
      }
    }

    // Scale factors for folding and faulting.
    float so = 1.0f/no;
    float sa = 1.0f/na;

    // Folding transforms
    Random r = new Random(2);
    int ng = 31;
    float[] g2 = mul(n2-1,randfloat(r,ng));
    float[] g3 = mul(n3-1,randfloat(r,ng));
    float[] sg = clip(4.0f,8.0f,mul(0.5f*(max(n2,n3)-1),randfloat(r,ng)));
    float[] hg = mul(neg(sg),randfloat(r,ng));
    T1 s1 = new Linear1(0.0f,so*2.0f/n1);
    T2 s2 = new Gaussians2(g2,g3,sg,hg);
    VerticalShear3 shear = new VerticalShear3(s1,s2);

    // Faulting transforms
    float r1a = 0.0f*n1, r2a = 0.4f*n2, r3a = 0.5f*n3;
    float r1b = 0.0f*n1, r2b = 0.1f*n2, r3b = 0.3f*n3;
    float r1c = 0.3f*n1, r2c = 0.7f*n2, r3c = 0.5f*n3;
    float r1d = 0.1f*n1, r2d = 0.5f*n2, r3d = 0.5f*n3;
    float phia =  10.0f, thetaa = 75.0f; if (conjugate) phia += 180.0f;
    float phib =  10.0f, thetab = 75.0f;
    float phic = 190.0f, thetac = 75.0f;
    float thetad = 75.0f;
    float[] c1 = {0.0f}, c2 = {0.0f}, sc = {20.0f}, hc = {sa*5.0f};
    T2 throwa = new Linear2(0.0f,sa*0.1f,0.0f,0.0f,0.0f,0.0f);
    T2 throwb = new Linear2(0.0f,sa*0.1f,0.0f,0.0f,0.0f,0.0f);
    T2 throwc = new Gaussians2(c1,c2,sc,hc);
    T1 throwd = new Linear1(0.0f,sa*0.1f);
    PlanarFault3 faulta = new PlanarFault3(r1a,r2a,r3a,phia,thetaa,throwa);
    PlanarFault3 faultb = new PlanarFault3(r1b,r2b,r3b,phib,thetab,throwb);
    PlanarFault3 faultc = new PlanarFault3(r1c,r2c,r3c,phic,thetac,throwc);
    ConicalFault3 faultd = new ConicalFault3(r1d,r2d,r3d,thetad,throwd);

    // Reflectivity or impedance.
    float[][][][] p = makeReflectivityWithNormals(m1,n2,n3);
    p = addChannels(p);
    if (impedance)
      p = impedanceFromReflectivity(p);

    // Apply the deformation sequence.
    for (int js=0; js<ns; ++js) {
      if (sequence.charAt(js)=='O') {
        p = apply(shear,p);
      } else if (sequence.charAt(js)=='A') {
        if (nplanar>0) p = apply(faulta,p);
        if (nplanar>1) p = apply(faultb,p);
        if (nplanar>2) p = apply(faultc,p);
        if (conical) p = apply(faultd,p);
      }
    }

    // Wavelet and noise.
    if (wavelet)
      p = addWavelet(0.15,p);
    p[0] = mul(1.0f/rms(p[0]),p[0]);
    p[0] = addNoise(noise,p[0]);

    // Slopes.
    p[2] = neg(div(p[2],p[1]));
    p[3] = neg(div(p[3],p[1]));

    // Trim the image.
    p[0] = copy(n1,n2,n3,p[0]);
    p[2] = copy(n1,n2,n3,p[2]);
    p[3] = copy(n1,n2,n3,p[3]);

    return new float[][][][]{p[0],p[2],p[3]};
  }

  /**
   * Returns a fake 3D seismic image with slopes.
   * @param sequence string of 'O' and 'A' for fOlding and fAulting.
   * @param nplanar number of planar faults.
   * @param conjugate true, to make one of the faults a conjugate fault.
   * @param conical true, to include a conical fault.
   * @param impedance true, for impedance instead of reflectivity.
   * @param wavelet true, for wavelet; false, for no wavelet.
   * @param noise rms of noise (relative to signal) added to the image.
   * @return array of arrays {f,p2,p3} with image f and slopes p2 and p3.
   */
  public static float[][][][] seismicAndSlopes3d2015A(
      String sequence, int nplanar, boolean conjugate, boolean conical,
      boolean impedance, boolean wavelet, boolean lateralViriation, double noise) {
    int n1 = 121;
    int n2 = 152;
    int n3 = 153;
    int m1 = n1+50;

    // Number of episodes in deformation sequence of folding and faulting.
    int ns = sequence.length();
    int no = 0;
    int na = 0;
    for (int js=0; js<ns; ++js) {
      if (sequence.charAt(js)=='O') {
        ++no;
      } else if (sequence.charAt(js)=='A') {
        ++na;
      }
    }

    // Scale factors for folding and faulting.
    float so = 1.0f/no;
    float sa = 1.0f/na;

    // Folding transforms
    Random r = new Random(2);
    int ng = 31;
    float[] g2 = mul(n2-1,randfloat(r,ng));
    float[] g3 = mul(n3-1,randfloat(r,ng));
    float[] sg = clip(4.0f,8.0f,mul(0.5f*(max(n2,n3)-1),randfloat(r,ng)));
    float[] hg = mul(neg(sg),randfloat(r,ng));
    T1 s1 = new Linear1(0.0f,so*2.0f/n1);
    T2 s2 = new Gaussians2(g2,g3,sg,hg);
    VerticalShear3X shearX = new VerticalShear3X(s1,s2);

    // Faulting transforms
    float r1a = 0.0f*n1, r2a = 0.5f*n2, r3a = 0.5f*n3;
    float r1b = 0.0f*n1, r2b = 0.2f*n2, r3b = 0.3f*n3;
    float r1c = 0.3f*n1, r2c = 0.8f*n2, r3c = 0.5f*n3;
    float phia =  10.0f, thetaa = 75.0f; if (conjugate) phia += 180.0f;
    float phib =  10.0f, thetab = 75.0f;
    float phic = 190.0f, thetac = 75.0f;
    float[] c1 = {0.0f}, c2 = {0.0f}, sc = {20.0f}, hc = {sa*8.0f};
    T2 throwa = new Linear2(0.0f,sa*0.1f,0.0f,0.0f,0.0f,0.0f);
    T2 throwb = new Linear2(0.0f,sa*0.1f,0.0f,0.0f,0.0f,0.0f);
    T2 throwc = new Gaussians2(c1,c2,sc,hc);
    PlanarFault3 faulta = new PlanarFault3(r1a,r2a,r3a,phia,thetaa,throwa);
    PlanarFault3 faultb = new PlanarFault3(r1b,r2b,r3b,phib,thetab,throwb);
    PlanarFault3 faultc = new PlanarFault3(r1c,r2c,r3c,phic,thetac,throwc);

    // Reflectivity or impedance.
    float[][] hz = new float[n3][n2];
    float[][][][] vd = makeVelocityAndDensityX(m1, n2, n3, lateralViriation);
    addChannelsXX(vd,hz);
    float[][][] v1 = copy(vd[0]);
    float[][][] v2 = copy(vd[0]);
    float[][][] d1 = copy(vd[1]);
    float[][][] d2 = copy(vd[1]);
    float[][][][] p = makeReflectivityWithNormals(vd);
    float[][][][] q = makeReflectivityWithNormals(vd);
    if (impedance)
      p = impedanceFromReflectivity(p);

    // Apply the deformation sequence.
    for (int js=0; js<ns; ++js) {
      if (sequence.charAt(js)=='O') {
        //p = apply(shear,p);
        q = apply(shearX,q);
        p = combine(32,p,q);
        //v1 = apply(shear,v1);
        v2 = apply(shearX,v2);
        v1 = combine(32,v1,v2);
        //d1 = apply(shear,d1);
        d2 = apply(shearX,d2);
        d1 = combine(32,d1,d2);
        apply(n1,shearX,hz);
        resetVelocityDensity(lateralViriation,32,p[0],v1,d1);
      } else if (sequence.charAt(js)=='A') {
        if (nplanar>0) {
          p  = apply(faulta,p); 
          v1 = apply(faulta,v1);
          d1 = apply(faulta,d1);
          apply(n1,faulta,hz);
        }
        if (nplanar>1) {
          p  = apply(faultb,p);
          v1 = apply(faultb,v1);
          d1 = apply(faultb,d1);
          apply(n1,faultb,hz);
        }
        if (nplanar>2) {
          p  = apply(faultc,p);
          v1 = apply(faultc,v1);
          d1 = apply(faultc,d1);
          apply(n1,faultc,hz);
        }
      }
    }

    float[][][] pc = new float[n3][n2][n1];
    float[][][] rc = new float[n3][n2][n1];
    float[][][] rn = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      float vi = v1[i3][i2][i1];
      float di = d1[i3][i2][i1];
      float vm = v1[i3][i2][i1-1];
      float dm = d1[i3][i2][i1-1];
      pc[i3][i2][i1] = vi*di;
      rc[i3][i2][i1] = 0.5f*(log(vi*di)-log(vm*dm));
    }}}
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float vi = v1[i3][i2][0];
      float di = d1[i3][i2][0];
      pc[i3][i2][0] = vi*di;
      rc[i3][i2][0] = rc[i3][i2][1];
    }}
    copy(rc,q[0]);
    copy(p[1],q[1]);
    copy(p[2],q[2]);
    copy(p[3],q[3]);
    rn = addNoise(noise,rc);
    // Wavelet and noise.
    if (wavelet) {
      q = addWavelet(0.15,q);
      p = addWavelet(0.15,p);
      p[0] = addNoise(noise,p[0]);
    }
    q[0] = mul(1.0f/rms(q[0]),q[0]);
    p[0] = mul(1.0f/rms(p[0]),p[0]);

    p[0] = copy(n1,n2,n3,p[0]);
    q[0] = copy(n1,n2,n3,q[0]);
    float[][][] hs = new float[][][]{hz};
    return new float[][][][]{q[0],p[0],rn,pc,hs};
  }

  /**
   * Returns a fake 3D seismic image with slopes.
   * @param sequence string of 'O' and 'A' for fOlding and fAulting.
   * @param nplanar number of planar faults.
   * @param conjugate true, to make one of the faults a conjugate fault.
   * @param conical true, to include a conical fault.
   * @param impedance true, for impedance instead of reflectivity.
   * @param wavelet true, for wavelet; false, for no wavelet.
   * @param noise rms of noise (relative to signal) added to the image.
   * @return array of arrays {f,p2,p3} with image f and slopes p2 and p3.
   */
  public static float[][][][] seismicAndSlopes3d2015B(
      String sequence, int nplanar, boolean conjugate, boolean conical,
      boolean impedance, boolean wavelet, boolean lateralViriation, double noise) {
    int n1 = 121;
    int n2 = 152;
    int n3 = 153;
    int m1 = n1+50;

    // Number of episodes in deformation sequence of folding and faulting.
    int ns = sequence.length();
    int no = 0;
    int na = 0;
    for (int js=0; js<ns; ++js) {
      if (sequence.charAt(js)=='O') {
        ++no;
      } else if (sequence.charAt(js)=='A') {
        ++na;
      }
    }

    // Scale factors for folding and faulting.
    float so = 1.0f/no;
    float sa = 1.0f/na;

    // Folding transforms
    Random r = new Random(2);
    int ng = 31;
    float[] g2 = mul(n2-1,randfloat(r,ng));
    float[] g3 = mul(n3-1,randfloat(r,ng));
    float[] sg = clip(4.0f,8.0f,mul(0.5f*(max(n2,n3)-1),randfloat(r,ng)));
    float[] hg = mul(neg(sg),randfloat(r,ng));
    T1 s1 = new Linear1(0.0f,so*2.0f/n1);
    T2 s2 = new Gaussians2(g2,g3,sg,hg);
    VerticalShear3X shearX = new VerticalShear3X(s1,s2);

    // Faulting transforms
    float r1a = 0.0f*n1, r2a = 0.5f*n2, r3a = 0.5f*n3;
    float r1b = 0.0f*n1, r2b = 0.2f*n2, r3b = 0.3f*n3;
    float r1c = 0.3f*n1, r2c = 0.8f*n2, r3c = 0.5f*n3;
    float phia =  10.0f, thetaa = 75.0f; if (conjugate) phia += 180.0f;
    float phib =  10.0f, thetab = 75.0f;
    float phic = 190.0f, thetac = 75.0f;
    float[] c1 = {0.0f}, c2 = {0.0f}, sc = {20.0f}, hc = {sa*8.0f};
    T2 throwa = new Linear2(0.0f,sa*0.1f,0.0f,0.0f,0.0f,0.0f);
    T2 throwb = new Linear2(0.0f,sa*0.1f,0.0f,0.0f,0.0f,0.0f);
    T2 throwc = new Gaussians2(c1,c2,sc,hc);
    PlanarFault3 faulta = new PlanarFault3(r1a,r2a,r3a,phia,thetaa,throwa);
    PlanarFault3 faultb = new PlanarFault3(r1b,r2b,r3b,phib,thetab,throwb);
    PlanarFault3 faultc = new PlanarFault3(r1c,r2c,r3c,phic,thetac,throwc);

    // Reflectivity or impedance.
    float[][][][] vd = makeVelocityAndDensityX(m1, n2, n3, lateralViriation);
    float[][] hz = new float[n3][n2];
    addChannelsXX(vd,hz);
    float[][][] v1 = copy(vd[0]);
    float[][][] v2 = copy(vd[0]);
    float[][][] d1 = copy(vd[1]);
    float[][][] d2 = copy(vd[1]);
    float[][][][] p = makeReflectivityWithNormals(vd);
    float[][][][] q = makeReflectivityWithNormals(vd);
    if (impedance)
      p = impedanceFromReflectivity(p);
    float[][][] pc = new float[n3][n2][n1];
    float[][][] rc = new float[n3][n2][n1];
    float[][][] rn = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      float vi = v1[i3][i2][i1];
      float di = d1[i3][i2][i1];
      float vm = v1[i3][i2][i1-1];
      float dm = d1[i3][i2][i1-1];
      pc[i3][i2][i1] = vi*di;
      rc[i3][i2][i1] = 0.5f*(log(vi*di)-log(vm*dm));
    }}}
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float vi = v1[i3][i2][0];
      float di = d1[i3][i2][0];
      pc[i3][i2][0] = vi*di;
      rc[i3][i2][0] = rc[i3][i2][1];
    }}
    float[][][] rct = copy(rc);
    float[][][] pct = copy(pc);
    // Apply the deformation sequence.
    for (int js=0; js<ns; ++js) {
      if (sequence.charAt(js)=='O') {
        q = apply(shearX,q);
        p = combine(32,p,q);
        rct = apply(shearX,rct);
        pct = apply(shearX,pct);
        rc  = combine(32,rc,rct);
        pc  = combine(32,pc,pct);
        /*
        v2 = apply(shearX,v2);
        v1 = combine(32,v1,v2);
        d2 = apply(shearX,d2);
        d1 = combine(32,d1,d2);
        */
        resetReflectivityImpedance(lateralViriation,32,p[0],rc,pc);
      } else if (sequence.charAt(js)=='A') {
        if (nplanar>0) {
          p  = apply(faulta,p); 
          rc = apply(faulta,rc);
          pc = apply(faulta,pc);
          //v1 = apply(faulta,v1);
          //d1 = apply(faulta,d1);
        }
        if (nplanar>1) {
          p  = apply(faultb,p);
          rc = apply(faultb,rc);
          pc = apply(faultb,pc);
          //v1 = apply(faultb,v1);
          //d1 = apply(faultb,d1);
        }
        if (nplanar>2) {
          p  = apply(faultc,p);
          rc = apply(faultc,rc);
          pc = apply(faultc,pc);
          //v1 = apply(faultc,v1);
          //d1 = apply(faultc,d1);
        }
      }
    }

    copy(rc,q[0]);
    copy(p[1],q[1]);
    copy(p[2],q[2]);
    copy(p[3],q[3]);
    rn = addNoise(noise,rc);
    // Wavelet and noise.
    if (wavelet) {
      q = addWavelet(0.15,q);
      p = addWavelet(0.15,p);
      p[0] = addNoise(noise,p[0]);
    }
    q[0] = mul(1.0f/rms(q[0]),q[0]);
    p[0] = mul(1.0f/rms(p[0]),p[0]);

    p[0] = copy(n1,n2,n3,p[0]);
    q[0] = copy(n1,n2,n3,q[0]);
    return new float[][][][]{q[0],p[0],rn,pc};
  }



  /**
   * Returns a fake noisy 2D seismic image with noise-free slopes.
   * The image and slopes are designed to be used to test methods for
   * estimating slopes of locally linear features apparent in 2D seismic
   * images. The image contains sinusoidal folding, horizontal and dipping
   * layers, two unconformities, and two intersecting faults with throws that
   * increase linearly with depth. The rms noise-to-signal ratio is 0.5.
   * @return array {f,p} of fake seismic image and corresponding slopes.
   */
  public static float[][][] seismicAndSlopes2d2014A() {
    return seismicAndSlopes2d2014A(0.5);
  }

  /**
   * Returns a fake 2D density and velocity models.
   * @param noise rms of noise (relative to signal) added to the image.
   */
  public static float[][][] densityAndVelocity2d(double noise, float[][][] rv) {
    int n1 = 501;
    int n2 = 501;
    enhanceContrast(rv);
    float[][][] p = densityAndVelocityFromLogs(n1,n2,rv);
    float[][][] q = densityAndVelocityFromLogs(n1,n2,rv);
    Sinusoidal2 fold = new Sinusoidal2(0.0f,0.04f,2.0e-2f,5.0e-5f);
    VerticalShear2 shear = new VerticalShear2(new Linear1(0.0f,0.1f));
    float[][] rh = apply(fold,p[0]);
    rh = combine(n1/3,q[0],rh);
    rh = apply(shear,rh);
    rh = combine(n1/6,q[0],rh);

    float[][] vp = apply(fold,p[1]);
    vp = combine(n1/3,q[1],vp);
    vp = apply(shear,vp);
    vp = combine(n1/6,q[1],vp);

    float[][] vs = apply(fold,p[2]);
    vs = combine(n1/3,q[2],vs);
    vs = apply(shear,vs);
    vs = combine(n1/6,q[2],vs);
    return new float[][][]{rh,vp,vs};
  }

  public static float[][][] densityAndVelocity2dx(double noise, float[][][] rv) {
    int n1 = 273;
    int n2 = 501;
    //enhanceContrast(rv);
    float[][][] p = densityAndVelocityFromLogs(n1,n2,rv);
    float[][][] q = densityAndVelocityFromLogs(n1,n2,rv);
    Sinusoidal2 fold = new Sinusoidal2(0.0f,0.04f,2.0e-2f,5.0e-5f);
    VerticalShear2 shear = new VerticalShear2(new Linear1(0.0f,0.1f));
    float[][] rh = apply(fold,p[0]);
    rh = combine(n1/3,q[0],rh);
    rh = apply(shear,rh);
    rh = combine(n1/6,q[0],rh);

    float[][] vp = apply(fold,p[1]);
    vp = combine(n1/3,q[1],vp);
    vp = apply(shear,vp);
    vp = combine(n1/6,q[1],vp);

    float[][] vs = apply(fold,p[2]);
    vs = combine(n1/3,q[2],vs);
    vs = apply(shear,vs);
    vs = combine(n1/6,q[2],vs);
    return new float[][][]{rh,vp,vs};
  }


  private static void enhanceContrast(float[][][] rv) {
    int n3 = rv.length;
    int n2 = rv[0].length;
    int n1 = rv[0][0].length;
    float[][][] rg = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply1XX(rv,rg);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      rv[i3][i2][i1] = rv[i3][i2][i1-1]+2f*rg[i3][i2][i1];
    }}}
  }


  /**
   * Returns a fake 2D seismic image with slopes.
   * The fake seismic image contains sinusoidal folding, horizontal and
   * dipping layers, two unconformities, and two intersecting faults with
   * throws that increase linearly with depth. While the image may have
   * a specified amount of additive noise, the slopes are noise-free.
   * @param noise rms of noise (relative to signal) added to the image.
   */
  public static float[][][] seismicAndSlopes2d2014A(double noise) {
    int n1 = 501;
    int n2 = 501;
    float[][][] p = makeReflectivityWithNormals(n1,n2);
    float[][][] q = makeReflectivityWithNormals(n1,n2);
    float[][][] r = makeReflectivityWithNormals(n1,n2);
    Linear1 throw1 = new Linear1(0.0f,0.10f);
    Linear1 throw2 = new Linear1(0.0f,0.10f);
    LinearFault2 fault1 = new LinearFault2(0.0f,n2*0.2f, 15.0f,throw1);
    LinearFault2 fault2 = new LinearFault2(0.0f,n2*0.4f,-15.0f,throw2);
    Sinusoidal2 fold = new Sinusoidal2(0.0f,0.05f,1.0e-4f,2.0e-4f);
    VerticalShear2 shear = new VerticalShear2(new Linear1(0.0f,0.05f));
    p = apply(fold,p);
    p = combine(n1/3,q,p);
    p = apply(shear,p);
    p = combine(n1/6,r,p);
    p = apply(fault1,p);
    p = apply(fault2,p);
    p = addWavelet(0.1,p);
    p[0] = addNoise(noise,p[0]);
    p[1] = neg(div(p[2],p[1]));
    return new float[][][]{p[0],p[1]};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static SincInterpolator _si = new SincInterpolator();
  static {
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
  }

  private static float[][] addNoise(double nrms, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    Random r = new Random(1);
    float[][] g = mul(2.0f,sub(randfloat(r,n1,n2),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2.0);
    rgf.apply10(g,g); // 1st derivative enhances high-frequencies
    g = mul(g,(float)nrms*rms(f)/rms(g));
    return add(f,g);
  }

  private static float[][][] addNoise(double nrms, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    Random r = new Random(1); // 31415
    float[][][] g = mul(2.0f,sub(randfloat(r,n1,n2,n3),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply100(g,g); // 1st derivative enhances high-frequencies
    g = mul(g,(float)nrms*rms(f)/rms(g));
    return add(f,g);
  }

  /**
   * Coordinates, Jacobians, and transforms.
   */
  private static class C2 {
    float c1,c2;
    C2(float c1, float c2) { 
      this.c1 = c1; 
      this.c2 = c2; 
    }
  }
  private static class C3 {
    float c1,c2,c3;
    C3(float c1, float c2, float c3) {
      this.c1 = c1;
      this.c2 = c2;
      this.c3 = c3;
    }
  }
  private static class D2 {
    float d11,d12,
          d21,d22;
    D2(float d11, float d12,
       float d21, float d22) { 
      this.d11 = d11;  this.d12 = d12; 
      this.d21 = d21;  this.d22 = d22; 
    }
  }
  private static class D3 {
    float d11,d12,d13,
          d21,d22,d23,
          d31,d32,d33;
    D3(float d11, float d12, float d13,
       float d21, float d22, float d23,
       float d31, float d32, float d33) { 
      this.d11 = d11;  this.d12 = d12;  this.d13 = d13;
      this.d21 = d21;  this.d22 = d22;  this.d23 = d23;
      this.d31 = d31;  this.d32 = d32;  this.d33 = d33;
    }
  }
  private interface T1 {
    float f(float x);
    float df(float x);
  }
  private interface T2 {
    C2 f(float x1, float x2);
    D2 df(float x1, float x2);
  }
  private interface T3 {
    C3 f(float x1, float x2, float x3);
    D3 df(float x1, float x2, float x3);
  }

  /**
   * Applies a 2D coordinate transform to an image.
   * @param t coordinate transform f(x).
   * @param p input image p(x).
   * @return transformed image q(x) = p(f(x)).
   */
  private static float[][] apply(T2 t, float[][] p) {
    int n1 = p[0].length;
    int n2 = p.length;
    float[][] q = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        C2 f = t.f(i1,i2);
        float f1 = f.c1;
        float f2 = f.c2;
        q[i2][i1] = _si.interpolate(n1,1.0,0.0,n2,1.0,0.0,p,f1,f2);
      }
    }
    return q;
  }

  private static void apply(
    int n1, T3 t, float[][] hz) {
    int n3 = hz.length;
    int n2 = hz[0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float[] x = new float[n1];
      float[] y = new float[n1];
      for (int i1=0; i1<n1; ++i1) {
        C3 f = t.f(i1,i2,i3);
        y[i1] = i1;
        x[i1] = f.c1;
      }
      CubicInterpolator ci = new CubicInterpolator(x,y);
      hz[i3][i2] = ci.interpolate(hz[i3][i2]);
    }}
  }


  /**
   * Applies a 2D coordinate transform to an image and normal vectors.
   * The input array {p0,p1,p2} contains the input image p0 and 1st and 2nd
   * components of normal vectors, p1 and p2. The returned array {q0,q1,q2}
   * contains the corresponding transformed image and normal vectors. All
   * normal vectors are unit vectors.
   * @param t coordinate transform f(x).
   * @param p input image p(x) and normal vectors.
   * @return transformed image q(x) = p(f(x) and normal vectors.
   */
  private static float[][][] apply(T2 t, float[][][] p) {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][] q0 = apply(t,p[0]);
    float[][] q1 = apply(t,p[1]);
    float[][] q2 = apply(t,p[2]);
    float[][] q3 = apply(t,p[3]);
    float[][] q4 = apply(t,p[4]);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        D2 d = t.df(i1,i2);
        float q3i = d.d11*q3[i2][i1]+d.d21*q4[i2][i1];
        float q4i = d.d12*q3[i2][i1]+d.d22*q4[i2][i1];
        float qsi = 1.0f/sqrt(q3i*q3i+q4i*q4i);
        q3[i2][i1] = q3i*qsi;
        q4[i2][i1] = q4i*qsi;
      }
    }
    return new float[][][]{q0,q1,q2,q3,q4};
  }

  /**
   * Applies a 3D coordinate transform to an image.
   * @param t coordinate transform f(x).
   * @param p input image p(x).
   * @return transformed image q(x) = p(f(x)).
   */
  private static float[][][] apply(final T3 t, final float[][][] p) {
    final int n1 = p[0][0].length;
    final int n2 = p[0].length;
    final int n3 = p.length;
    final float[][][] q = new float[n3][n2][n1];
    //for (int i3=0; i3<n3; ++i3) {
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          C3 f = t.f(i1,i2,i3);
          float f1 = f.c1;
          float f2 = f.c2;
          float f3 = f.c3;
          q[i3][i2][i1] = _si.interpolate(
              n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,p,f1,f2,f3);
        }
      }
    }});
    //}
    return q;
  }

  /**
   * Applies a 3D coordinate transform to an image and normal vectors.
   * The input array {p0,p1,p2,p3} contains the input image p0 and 1st, 2nd,
   * and 3rd components of normal vectors, p1, p2, and p3. The returned array
   * {q0,q1,q2,q3} contains the corresponding transformed image and normal
   * vectors. All normal vectors are unit vectors.
   * @param t coordinate transform f(x).
   * @param p input image p(x) and normal vectors.
   * @return transformed image q(x) = p(f(x) and normal vectors.
   */
  private static float[][][][] apply(final T3 t, final float[][][][] p) {
    final int n1 = p[0][0][0].length;
    final int n2 = p[0][0].length;
    final int n3 = p[0].length;
    final float[][][] q0 = apply(t,p[0]);
    final float[][][] q1 = apply(t,p[1]);
    final float[][][] q2 = apply(t,p[2]);
    final float[][][] q3 = apply(t,p[3]);
    //for (int i3=0; i3<n3; ++i3) {
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          D3 d = t.df(i1,i2,i3);
          float t1i = q1[i3][i2][i1];
          float t2i = q2[i3][i2][i1];
          float t3i = q3[i3][i2][i1];
          float q1i = d.d11*t1i+d.d21*t2i+d.d31*t3i;
          float q2i = d.d12*t1i+d.d22*t2i+d.d32*t3i;
          float q3i = d.d13*t1i+d.d23*t2i+d.d33*t3i;
          float qsi = 1.0f/sqrt(q1i*q1i+q2i*q2i+q3i*q3i);
          q1[i3][i2][i1] = q1i*qsi;
          q2[i3][i2][i1] = q2i*qsi;
          q3[i3][i2][i1] = q3i*qsi;
        }
      }
    }});
    //}
    return new float[][][][]{q0,q1,q2,q3};
  }

  /**
   * A linear 1D coordinate mapping f(x) = a0+a1*x.
   */
  private static class Linear1 implements T1 {
    public Linear1(float a0, float a1) { 
      _a0 = a0; _a1 = a1;
    }
    public float f(float x) { 
      return _a0+_a1*x; 
    }
    public float df(float x) { 
      return _a1; 
    }
    private float _a0,_a1;
  }

  /**
   * A linear 2D coordinate mapping f(x) = a0+a1*x1+a2*x2.
   * In this mapping, f, x, a0, a1, and a2 are vectors, each with two
   * components.
   */
  private static class Linear2 implements T2 {
    public Linear2(
        float a10, float a11, float a12,
        float a20, float a21, float a22) {
      _a10 = a10; _a11 = a11; _a12 = a12;
      _a20 = a20; _a21 = a21; _a22 = a22;
    }
    public C2 f(float x1, float x2) {
      float f1 = _a10+_a11*x1+_a12*x2;
      float f2 = _a20+_a21*x1+_a22*x2;
      return new C2(f1,f2);
    }
    public D2 df(float x1, float x2) {
      return new D2(_a11,_a12,_a21,_a22);
    }
    private float _a10,_a11,_a12,_a20,_a21,_a22;
  }

  /**
   * A linear fault in a 2D image.
   */
  private static class LinearFault2 implements T2 {

    /**
     * Constructs a linear fault.
     * @param fx1 coordinate x1 of a reference point on the fault.
     * @param fx2 coordinate x2 of a reference point on the fault.
     * @param ftheta fault dip, measured in degrees from vertical.
     * @param fthrow fault throw, a function of coordinate x1.
     */
    public LinearFault2(float fx1, float fx2, float ftheta, T1 fthrow) {

      // Reference point (the origin in fault-line coordinates).
      _r1 = fx1;
      _r2 = fx2;

      // Tangent of fault dip.
      float rtheta = toRadians(ftheta);
      _ttheta = tan(rtheta);

      // Fault normal vector.
      float ctheta = cos(rtheta);
      float stheta = sin(rtheta);
      _u1 = -stheta;
      _u2 =  ctheta;

      // Ensure vertical component of normal vector is non-positive.
      if (_u1>0.0f) {
        _u1 = -_u1;
        _u2 = -_u2;
      }

      // Constant needed to locate points with respect to plane.
      _u0 = -(fx1*_u1+fx2*_u2);

      // Fault throw.
      _t1 = fthrow;
    }
    public C2 f(float x1, float x2) {
      if (faulted(x1,x2)) {
        x1 -= _r1;
        x2 -= _r2;
        float t = _t1.f(x1);
        x1 -= t;
        x2 -= t*_ttheta;
        x1 += _r1;
        x2 += _r2;
      }
      return new C2(x1,x2);
    }
    public D2 df(float x1, float x2) {
      float d11 = 1.0f, d12 = 0.0f,
            d21 = 0.0f, d22 = 1.0f;
      if (faulted(x1,x2)) {
        x1 -= _r1;
        x2 -= _r2;
        float dt = _t1.df(x1);
        d11 -= dt;
        d21 -= dt*_ttheta;
      }
      return new D2(d11,d12,
                    d21,d22);
    }
    private float _r1,_r2;
    private float _ttheta;
    private float _u0,_u1,_u2;
    private T1 _t1;
    private boolean faulted(float x1, float x2) {
      return _u0+_u1*x1+_u2*x2>=0.0f;
    }
  }

  /**
   * A planar fault in a 3D image.
   */
  private static class PlanarFault3 implements T3 {

    /**
     * Constructs a planar fault. Fault strike and dip angles define
     * a vector normal to the fault. With this fault normal vector and
     * a reference point that lies on the fault, we can locate points
     * with respect to the fault plane. Points on one side of the fault
     * will be displaced relative to points on the other side. The
     * reference point corresponds to the origin of the 2D coordinate
     * system in which fault throws are specified.
     * @param fx1 coordinate x1 of a reference point on the fault.
     * @param fx2 coordinate x2 of a reference point on the fault.
     * @param fx3 coordinate x3 of a reference point on the fault.
     * @param fphi fault strike, measured in degrees from x3 axis.
     * @param ftheta fault dip, measured in degrees from horizontal.
     * @param fthrow fault throw, function of fault-plane coordinates.
     */
    public PlanarFault3(
        float fx1, float fx2, float fx3,
        float fphi, float ftheta, T2 fthrow) {

      // Reference point (the origin in fault-plane coordinates).
      _r1 = fx1;
      _r2 = fx2;
      _r3 = fx3;

      // Cosine and sine of fault strike.
      float rphi = toRadians(fphi);
      _cphi = cos(rphi);
      _sphi = sin(rphi);

      // Tangent of fault dip.
      float rtheta = toRadians(ftheta);
      _ottheta = 1.0f/tan(rtheta);

      // Fault normal vector.
      float[] w = faultNormalVectorFromStrikeAndDip(fphi,ftheta);
      _w1 = w[0];
      _w2 = w[1];
      _w3 = w[2];

      // Constant needed to locate points with respect to plane.
      _w0 = -(fx1*_w1+fx2*_w2+fx3*_w3);

      // Rotation matrix used to align fault strike with axis 3.
      _r22 =  _cphi;
      _r23 = -_sphi;
      _r32 =  _sphi;
      _r33 =  _cphi;

      // Fault throw.
      _t2 = fthrow;
    }
    public C3 f(float x1, float x2, float x3) {
      if (faulted(x1,x2,x3)) {
        x1 -= _r1;
        x2 -= _r2;
        x3 -= _r3;
        float y1 = x1;
        float y2 = _r22*x2+_r23*x3;
        float y3 = _r32*x2+_r33*x3;
        C2 t = _t2.f(y1,y3);
        y1 -= t.c1;
        y2 -= t.c1*_ottheta;
        x1 = y1;
        x2 = _r22*y2+_r32*y3;
        x3 = _r23*y2+_r33*y3;
        x1 += _r1;
        x2 += _r2;
        x3 += _r3;
      }
      return new C3(x1,x2,x3);
    }
    public D3 df(float x1, float x2, float x3) {
      float d11 = 1.0f, d12 = 0.0f, d13 = 0.0f,
            d21 = 0.0f, d22 = 1.0f, d23 = 0.0f,
            d31 = 0.0f, d32 = 0.0f, d33 = 1.0f;
      if (faulted(x1,x2,x3)) {
        x1 -= _r1;
        x2 -= _r2;
        x3 -= _r3;
        float y1 = x1;
        float y3 = _r32*x2+_r33*x3;
        D2 dt = _t2.df(y1,y3);
        float y11 = -dt.d11;
        float y13 = -dt.d12;
        float y21 = y11*_ottheta;
        float y23 = y13*_ottheta;
        d11 += y11;        d12 += y13*_sphi;        d13 += y13*_cphi;
        d21 += y21*_cphi;  d22 += y23*_cphi*_sphi;  d23 += y23*_cphi*_cphi;
        d31 -= y21*_sphi;  d32 -= y23*_sphi*_sphi;  d33 -= y23*_cphi*_sphi;
      }
      return new D3(d11,d12,d13,
                    d21,d22,d23,
                    d31,d32,d33);
    }
    private float _r1,_r2,_r3;
    private float _cphi,_sphi,_ottheta;
    private float _w0,_w1,_w2,_w3;
    private float _r22,_r23,_r32,_r33;
    private T2 _t2;
    private boolean faulted(float x1, float x2, float x3) {
      return _w0+_w1*x1+_w2*x2+_w3*x3>=0.0f;
    }
  }

  /**
   * A conical fault in a 3D image.
   */
  private static class ConicalFault3 implements T3 {

    /**
     * Constructs a conical fault. The cone geometry is defined by the
     * location of its apex and by its dip angle. Throw is a function of depth
     * relative to the apex, where throw must be zero.
     * @param fa1 coordinate x1 of fault apex.
     * @param fa2 coordinate x2 of a reference point on the fault.
     * @param fa3 coordinate x3 of a reference point on the fault.
     * @param ftheta fault dip, measured in degrees from x1 axis.
     * @param fthrow fault throw.
     */
    public ConicalFault3(
        float fa1, float fa2, float fa3,
        float ftheta, T1 fthrow) {

      // Apex of cone, the origin in fault coordinates.
      _a1 = fa1;
      _a2 = fa2;
      _a3 = fa3;

      // Inverse tangent of fault dip.
      float rtheta = toRadians(ftheta);
      _ottheta = 1.0f/tan(rtheta);

      // Fault throw.
      _t1 = fthrow;
    }
    public C3 f(float x1, float x2, float x3) {
      if (faulted(x1,x2,x3)) {
        x1 -= _a1;
        x2 -= _a2;
        x3 -= _a3;
        float t = _t1.f(x1);
        float h = t*_ottheta;
        float xs = 1.0f/sqrt(x2*x2+x3*x3);
        x1 -= t;
        x2 -= h*x2*xs;
        x3 -= h*x3*xs;
        x1 += _a1;
        x2 += _a2;
        x3 += _a3;
      }
      return new C3(x1,x2,x3);
    }
    public D3 df(float x1, float x2, float x3) {
      float d11 = 1.0f, d12 = 0.0f, d13 = 0.0f,
            d21 = 0.0f, d22 = 1.0f, d23 = 0.0f,
            d31 = 0.0f, d32 = 0.0f, d33 = 1.0f;
      if (faulted(x1,x2,x3)) {
        x1 -= _a1;
        x2 -= _a2;
        x3 -= _a3;
        float t = _t1.f(x1);
        float dt = _t1.df(x1);
        float xs = 1.0f/sqrt(x2*x2+x3*x3);
        float xs3 = xs*xs*xs;
        d11 -= dt;
        d21 -= dt*_ottheta*x2*xs;
        d22 -= t*_ottheta*x3*x3*xs3;
        d23 += t*_ottheta*x2*x3*xs3;
        d31 -= dt*_ottheta*x3*xs;
        d32 += t*_ottheta*x2*x3*xs3;
        d33 -= t*_ottheta*x3*x3*xs3;
      }
      return new D3(d11,d12,d13,
                    d21,d22,d23,
                    d31,d32,d33);
    }
    private float _a1,_a2,_a3;
    private float _ottheta;
    private T1 _t1;
    private boolean faulted(float x1, float x2, float x3) {
      if (x1<=_a1)
        return false;
      float d1 = x1-_a1;
      float d2 = x2-_a2;
      float d3 = x3-_a3;
      return d2*d2+d3*d3>d1*d1*_ottheta*_ottheta;
    }
  }

  /**
   * Vertical shear of a 2D image.
   */
  private static class VerticalShear2 implements T2 {
    public VerticalShear2(T1 s1) {
      _s1 = s1;
    }
    public C2 f(float x1, float x2) {
      x1 -= _s1.f(x2);
      return new C2(x1,x2);
    }
    public D2 df(float x1, float x2) {
      float d12 = -_s1.df(x2);
      return new D2(1.0f,  d12,
                    0.0f, 1.0f);
    }
    private T1 _s1;
  }

  /**
   * Vertical shear of a 3D image.
   */
  private static class VerticalShear3 implements T3 {
    public VerticalShear3(T1 s1, T2 s2) {
      _s1 = s1;
      _s2 = s2;
    }
    public C3 f(float x1, float x2, float x3) {
      float c1 = _s1.f(x1);
      C2 c2 = _s2.f(x2,x3);
      x1 -= c1*c2.c1;
      return new C3(x1,x2,x3);
    }
    public D3 df(float x1, float x2, float x3) {
      float c1 = _s1.f(x1);
      float d1 = _s1.df(x1);
      C2 c2 = _s2.f(x2,x3);
      D2 d2 = _s2.df(x2,x3);
      float d11 = 1.0f-d1*c2.c1;
      float d12 = -c1*d2.d11;
      float d13 = -c1*d2.d12;
      return new D3( d11,  d12,  d13,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f);
    }
    private T1 _s1;
    private T2 _s2;
  }

  private static class VerticalShear3X implements T3 {
    public VerticalShear3X(T1 s1, T2 s2) {
      _s1 = s1;
      _s2 = s2;
    }
    public C3 f(float x1, float x2, float x3) {
      float c1 = _s1.f(x1);
      C2 c2 = _s2.f(x2,x3);
      x1 -= (c1*c2.c1+x2*0.15f);
      return new C3(x1,x2,x3);
    }
    public D3 df(float x1, float x2, float x3) {
      float c1 = _s1.f(x1);
      float d1 = _s1.df(x1);
      C2 c2 = _s2.f(x2,x3);
      D2 d2 = _s2.df(x2,x3);
      float d11 = 1.0f-d1*c2.c1;
      float d12 = -c1*d2.d11;
      float d13 = -c1*d2.d12;
      return new D3( d11,  d12,  d13,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f);
    }
    private T1 _s1;
    private T2 _s2;
  }


  /**
   * A sum of 2D Gaussians. Each Gaussian is centered at (c1,c2), with radius
   * sg, and a height of hg.
   */
  private static class Gaussians2 implements T2 {
    public Gaussians2(float[] g1, float[] g2, float[] sg, float[] hg) {
      _ng = g1.length;
      _g1 = copy(g1);
      _g2 = copy(g2);
      _ag = copy(hg);
      _bg = div(-0.5f,mul(sg,sg));
    }
    public C2 f(float x1, float x2) {
      float f1 = 0.0f;
      float f2 = 0.0f;
      for (int ig=0; ig<_ng; ++ig) {
        float d1 = x1-_g1[ig];
        float d2 = x2-_g2[ig];
        f1 += _ag[ig]*exp(_bg[ig]*(d1*d1+d2*d2));
      }
      return new C2(f1,f2);
    }
    public D2 df(float x1, float x2) {
      float d11 = 0.0f, d12 = 0.0f;
      float d21 = 0.0f, d22 = 0.0f;
      for (int ig=0; ig<_ng; ++ig) {
        float d1 = x1-_g1[ig];
        float d2 = x2-_g2[ig];
        float si = 2.0f*_ag[ig]*_bg[ig];
        d11 += si*d1*exp(_bg[ig]*(d1*d1+d2*d2));
        d12 += si*d2*exp(_bg[ig]*(d1*d1+d2*d2));
      }
      return new D2(d11,d12,d21,d22);
    }
    private int _ng;
    private float[] _g1,_g2,_ag,_bg;
  }

  /**
   * Sinusoidal folding in a 2D image.
   */
  private static class Sinusoidal2 implements T2 {
    public Sinusoidal2(float a1, float b1, float a2, float b2) {
      _a1 = a1;
      _b1 = b1;
      _a2 = a2;
      _b2 = b2;
    }
    public C2 f(float x1, float x2) {
      float f1 = x1-(_a1+_b1*x1)*sin((_a2+_b2*x2)*x2);
      float f2 = x2;
      return new C2(f1,f2);
    }
    public D2 df(float x1, float x2) {
      float d11 = 1.0f-_b1*sin((_a2+_b2*x2)*x2);
      float d12 = -(_a1+_b1*x1)*(_a2+2.0f*_b2*x2)*cos((_a2+_b2*x2)*x2);
      float d21 = 0.0f;
      float d22 = 1.0f;
      return new D2(d11,d12,d21,d22);
    }
    private float _a1,_b1,_a2,_b2;
  }

  private static float[][][] densityAndVelocityFromLogs(
    int n1, int n2, float[][][] rv) {
    float[] x = new float[4];
    x[0] = 0;
    x[1] = 180;
    x[2] = 360;
    x[3] = 500;
    float[][] rh = new float[n2][n1];
    float[][] vp = new float[n2][n1];
    float[][] vs = new float[n2][n1];
    for (int i1=0; i1<n1; ++i1) {
      float[] rhi = new float[4];
      float[] vpi = new float[4];
      float[] vsi = new float[4];
      rhi[0] = rv[0][0][i1];
      rhi[1] = rv[0][1][i1];
      rhi[2] = rv[0][2][i1];
      rhi[3] = rv[0][3][i1];

      vpi[0] = rv[1][0][i1];
      vpi[1] = rv[1][1][i1];
      vpi[2] = rv[1][2][i1];
      vpi[3] = rv[1][3][i1];

      vsi[0] = rv[2][0][i1];
      vsi[1] = rv[2][1][i1];
      vsi[2] = rv[2][2][i1];
      vsi[3] = rv[2][3][i1];
      CubicInterpolator cir = new CubicInterpolator(x,rhi);
      CubicInterpolator cip = new CubicInterpolator(x,vpi);
      CubicInterpolator cis = new CubicInterpolator(x,vsi);
      for (int i2=0; i2<n2; ++i2) {
        rh[i2][i1] = cir.interpolate(i2);
        vp[i2][i1] = cip.interpolate(i2);
        vs[i2][i1] = cis.interpolate(i2);
      }
    }
    float[][][] p = new float[5][n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      copy(rh[i2],p[0][i2]);
      copy(vp[i2],p[1][i2]);
      copy(vs[i2],p[2][i2]);
    }
    p[3] = fillfloat(1.0f,n1,n2);
    p[4] = fillfloat(0.0f,n1,n2);
    return p;
  }

  private static float[][][] makeReflectivityWithNormals(int n1, int n2) {
    Random random = new Random(31);
    float[] r = pow(mul(2.0f,sub(randfloat(random,n1),0.5f)),5.0f);
    float[][][] p = new float[3][n2][n1];
    for (int i2=0; i2<n2; ++i2)
      copy(r,p[0][i2]);
    p[1] = fillfloat(1.0f,n1,n2);
    p[2] = fillfloat(0.0f,n1,n2);
    return p;
  }

  private static float[][][][] makeReflectivityWithNormals(
      int n1, int n2, int n3) {
    Random random = new Random(31);
    float[] r = mul(0.1f,pow(mul(2.0f,sub(randfloat(random,n1),0.5f)),5.0f));
    float[][][][] p = new float[4][n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        copy(r,p[0][i3][i2]);
      }
    }
    p[1] = fillfloat(1.0f,n1,n2,n3);
    p[2] = fillfloat(0.0f,n1,n2,n3);
    p[3] = fillfloat(0.0f,n1,n2,n3);
    return p;
  }

  private static float[][][][] impedanceFromReflectivity(float[][][][] r) {
    int n1 = r[0][0][0].length;
    int n2 = r[0][0].length;
    int n3 = r[0].length;
    float[][][] z = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        z[i3][i2][0] = 1.0f;
        for (int i1=1; i1<n1; ++i1) {
          float ri = r[0][i3][i2][i1];
          z[i3][i2][i1] = z[i3][i2][i1-1]*(1.0f+ri)/(1.0f-ri);
        }
      }
    }
    return new float[][][][]{z,r[1],r[2],r[3]};
  }

  private static float[][][][] makeVelocityAndDensityX(
    int n1, int n2, int n3, boolean lateralViriation) {
    int dn1 = 13;
    int dn2 = 22;
    int dn3 = 20;
    Random random = new Random(43);
    float[] r = randfloat(random,n1);
    float[][][] v = new float[n3][n2][n1];
    float[][][] d = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(i1<=dn1) {
        v[i3][i2][i1] = 0.5f*r[i1]+3.0f;
        d[i3][i2][i1] = 0.2f*r[i1]+2.2f;
      }
      if(i1>dn1&&i1<=dn1+dn2) {
        v[i3][i2][i1] = 0.5f*r[i1]+3.5f;
        d[i3][i2][i1] = 0.2f*r[i1]+2.4f;
      }
      if(i1>dn1+dn2&&i1<=dn1+dn2+dn3) {
        v[i3][i2][i1] = 0.6f*r[i1]+4.0f;
        d[i3][i2][i1] = 0.2f*r[i1]+2.6f;
      }
      if(i1>dn1+dn2+dn3) {
        v[i3][i2][i1] = 0.8f*r[i1]+4.5f;
        d[i3][i2][i1] = 0.2f*r[i1]+2.8f;
      }
      v[i3][i2][i1] *= 1000f;
    }}}
    if(lateralViriation) {
      makeLateralViriation(200f,v);
      makeLateralViriation(.07f,d);
    }
    return new float[][][][]{v,d};
  }

  private static void makeLateralViriation(float df, float[][][] f) {
    int d1 = 15;
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int i3=0; i3<n3; ++i3) {
    for (int i1=0; i1<n1; i1+=d1*2) {
      int nk = 7;
      double dk = 10;
      double fk = 20;
      for (int k1=i1; k1<i1+d1 && k1<n1; ++k1) {
        float fi = f[i3][0][k1];
        float[] fu = fillfloat(fi,7);
        fu[0] = fi-df;
        fu[6] = fi+df;
        for (int i2=0; i2<n2; ++i2) {
          f[i3][i2][k1] = si.interpolate(nk,dk,fk,fu,i2);
        }
        fk += 2.0;
      }
    }}
    for (int i3=0; i3<n3; ++i3) {
    for (int i1=d1; i1<n1; i1+=d1*2) {
      int nk = 7;
      double dk = 10;
      double fk = 50;
      for (int k1=i1; k1<i1+d1 && k1<n1; ++k1) {
        float fi = f[i3][0][k1];
        float[] fu = fillfloat(fi,7);
        fu[0] = fi-df;
        fu[6] = fi+df;
        for (int i2=0; i2<n2; ++i2) {
          f[i3][i2][k1] = si.interpolate(nk,dk,fk,fu,i2);
        }
        fk -= 2.0;
      }
    }}

  }

  private static float[][][][] makeVelocityAndDensity(int n1, int n2, int n3) {
    int dn1 = round((n1-50)/4f)-20;
    int dn2 = dn1+20;
    int dn3 = dn1+20;
    int dn4 = n1-dn1-dn2-dn3;
    Random random = new Random(31);
    float[] z = new float[n1];
    float[] r = mul(0.05f,pow(mul(2.0f,sub(randfloat(random,n1),0.5f)),5.0f));
    z[0] = 6000.0f;
    for (int i1=1; i1<n1; ++i1) {
      float ri = r[i1];
       z[i1] = z[i1-1]*(1.0f+ri)/(1.0f-ri);
    }
    Random random1 = new Random(round((31f*dn1)/n1));
    Random random2 = new Random(round((31f*dn2)/n1));
    Random random3 = new Random(round((31f*dn3)/n1));
    Random random4 = new Random(round((31f*dn4)/n1));
    float[] dv1 = mul(1000f,mul(0.50f,add(randfloat(random1,dn1),5.0f)));
    float[] dv2 = mul(1000f,mul(0.45f,add(randfloat(random2,dn2),5.0f)));
    float[] dv3 = mul(1000f,mul(0.55f,add(randfloat(random3,dn3),5.0f)));
    float[] dv4 = mul(1000f,mul(0.50f,add(randfloat(random4,dn4),6.0f)));

    float[] d1 = new float[n1];
    float[] v1 = new float[n1];
    copy(dn1,0,dv1,0,v1);
    copy(dn2,0,dv2,dn1,v1);
    copy(dn3,0,dv3,dn1+dn2,v1);
    copy(dn4,0,dv4,dn1+dn2+dn3,v1);
    for (int i1=0; i1<n1; ++i1) {
       float di = z[i1]/v1[i1];
       d1[i1] = di;
    }

    float[][][] d = new float[n3][n2][n1];
    float[][][] v = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      copy(v1,v[i3][i2]);
      copy(d1,d[i3][i2]);
    }}
    return new float[][][][] {v,d};
  }

  private static void resetVelocityDensity(boolean lateralViriation,
    int dn1, float[][][] p, float[][][] v, float[][][] d) 
  {

    Random random = new Random(43);
    int n3 = v.length;
    int n2 = v[0].length;
    int n1 = v[0][0].length;
    float[] r = randfloat(random,n1);
    float[][][] ds = new float[n3][n2][dn1];
    float[][][] vs = new float[n3][n2][dn1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<dn1; ++i1) {
      ds[i3][i2][i1] = 0.2f*r[i1]+2.2f;
      vs[i3][i2][i1] = (0.5f*r[i1]+3.0f)*1000f;
    }}}
    if(lateralViriation) {
      makeLateralViriation(200f,vs);
      makeLateralViriation(.07f,ds);
    }
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<dn1; ++i1) {
      d[i3][i2][i1] = ds[i3][i2][i1];
      v[i3][i2][i1] = vs[i3][i2][i1];
    }}}
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<dn1; ++i1) {
      float vi = v[i3][i2][i1];
      float di = d[i3][i2][i1];
      float vp = v[i3][i2][i1+1];
      float dp = d[i3][i2][i1+1];
      p[i3][i2][i1] = (vp*dp-vi*di)/(vp*dp+vi*di);
    }}}

  }

  private static void resetReflectivityImpedance(boolean lateralViriation,
    int dn1, float[][][] p, float[][][] rf, float[][][] ip) 
  {
    Random random = new Random(43);
    int n3 = rf.length;
    int n2 = rf[0].length;
    int n1 = rf[0][0].length;
    float[] r = randfloat(random,n1);
    float[][][] ds = new float[n3][n2][dn1];
    float[][][] vs = new float[n3][n2][dn1];
    float[][][] rs = new float[n3][n2][dn1];
    float[][][] ps = new float[n3][n2][dn1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<dn1; ++i1) {
      ds[i3][i2][i1] = 0.2f*r[i1]+2.2f;
      vs[i3][i2][i1] = (0.5f*r[i1]+3.0f)*1000f;
    }}}
    if(lateralViriation) {
      makeLateralViriation(200f,vs);
      makeLateralViriation(.07f,ds);
    }
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<dn1; ++i1) {
      float vi = vs[i3][i2][i1];
      float di = ds[i3][i2][i1];
      float vm = vs[i3][i2][i1-1];
      float dm = ds[i3][i2][i1-1];
      ps[i3][i2][i1] = vi*di;
      rs[i3][i2][i1] = 0.5f*(log(vi*di)-log(vm*dm));
    }}}
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float vi = vs[i3][i2][0];
      float di = ds[i3][i2][0];
      ps[i3][i2][0] = vi*di;
      rs[i3][i2][0] = rs[i3][i2][1];
    }}

    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<dn1; ++i1) {
      rf[i3][i2][i1] = rs[i3][i2][i1];
      ip[i3][i2][i1] = ps[i3][i2][i1];
      p[i3][i2][i1]  = ps[i3][i2][i1];
    }}}
  }



  private static float[][][][] makeReflectivityWithNormals(
      float[][][][] vd) {
    float[][][] v = vd[0];
    float[][][] d = vd[1];
    int n3 = v.length;
    int n2 = v[0].length;
    int n1 = v[0][0].length;
    float[][][][] p = new float[4][n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1-1; ++i1) {
      float vi = v[i3][i2][i1];
      float di = d[i3][i2][i1];
      float vp = v[i3][i2][i1+1];
      float dp = d[i3][i2][i1+1];
      p[0][i3][i2][i1] = (vp*dp-vi*di)/(vp*dp+vi*di);
    }}}
    p[1] = fillfloat(1.0f,n1,n2,n3);
    p[2] = fillfloat(0.0f,n1,n2,n3);
    p[3] = fillfloat(0.0f,n1,n2,n3);
    return p;
  }


  private static float[][][][] addChannels(float[][][][] r) {
    int n2 = r[0][0].length;
    int n3 = r[0].length;
    float[][][] s = copy(r[0]);
    int k1 = 75;
    float[] c2 = {0.4f*n2,0.6f*n2,0.5f*n2,0.7f*n2};
    float[] c3 = {0.1f*n3,0.5f*n3,0.7f*n3,0.9f*n3};
    CubicInterpolator.Method method = CubicInterpolator.Method.SPLINE;
    CubicInterpolator ci = new CubicInterpolator(method,c3,c2);
    for (int i3=0; i3<n3; ++i3) {
      float x2 = ci.interpolate(i3);
      int j2 = (int)(x2+0.5);
      for (int i2=j2-5; i2<j2+5; ++i2) {
        if (0<i2 && i2<n2) {
          float d2 = i2-j2;
          float a = 0.1f*exp(-0.125f*d2*d2);
          s[i3][i2][k1-1] -= 0.5f*a;
          s[i3][i2][k1  ] += 1.0f*a;
          s[i3][i2][k1+1] -= 0.5f*a;
        }
      }
    }
    return new float[][][][]{s,r[1],r[2],r[3]};
  }

  private static void addChannelsX(float[][][][] r) {
    int n3 = r[0].length;
    int n2 = r[0][0].length;
    int k1 = 70;
    float[] c2 = {0.4f*n2,0.6f*n2,0.5f*n2,0.7f*n2};
    float[] c3 = {0.1f*n3,0.5f*n3,0.7f*n3,0.9f*n3};
    CubicInterpolator.Method method = CubicInterpolator.Method.SPLINE;
    CubicInterpolator ci = new CubicInterpolator(method,c3,c2);
    for (int i3=0; i3<n3; ++i3) {
      float x2 = ci.interpolate(i3);
      int j2 = (int)(x2+0.5);
      for (int i2=j2-6; i2<j2+6; ++i2) {
        if (0<i2 && i2<n2) {
          float d2 = i2-j2;
          float a = 0.1f*exp(-0.125f*d2*d2);
          r[0][i3][i2][k1-1] += 4500f*a;
          r[0][i3][i2][k1  ] += 9000f*a;
          r[0][i3][i2][k1+1] += 4500f*a;

          r[1][i3][i2][k1-1] += 1.5f*a;
          r[1][i3][i2][k1  ] += 3.0f*a;
          r[1][i3][i2][k1+1] += 1.5f*a;
        }
      }
    }
  }

  private static void addChannelsXX(float[][][][] r, float[][] hz) {
    int n3 = r[0].length;
    int n2 = r[0][0].length;
    int k1 = 70;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      hz[i3][i2] = k1+1;
    }}
    float[] c2 = {0.1f*n2,0.5f*n2,0.7f*n2,0.9f*n2};
    float[] c3 = {0.4f*n3,0.6f*n3,0.5f*n3,0.7f*n3};
    CubicInterpolator.Method method = CubicInterpolator.Method.SPLINE;
    CubicInterpolator ci = new CubicInterpolator(method,c2,c3);
    for (int i2=0; i2<n2; ++i2) {
      float x3 = ci.interpolate(i2);
      int j3 = (int)(x3+0.5);
      for (int i3=j3-5; i3<j3+5; ++i3) {
        if (0<i3 && i3<n3) {
          float d2 = i3-j3;
          float a = 0.1f*exp(-0.125f*d2*d2);
          r[0][i3][i2][k1-1] += 9000f*a;
          r[0][i3][i2][k1  ] += 9000f*a;
          r[0][i3][i2][k1+1] += 9000f*a;
          //r[1][i3][i2][k1-1] += 2.0f*a;
          //r[1][i3][i2][k1  ] += 3.0f*a;
          //r[1][i3][i2][k1+1] += 2.0f*a;
        }
      }
    }
  }


  private static float[][] combine(
    float depth, float[][] pa, float[][] pb) 
  {
    int n2 = pa.length;
    int n1 = pa[0].length;
    float[][] pc = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          pc[i2][i1] = (i1<depth)?pa[i2][i1]:pb[i2][i1];
        }
      }
    return pc;
  }

  private static float[][][] combine(
    float depth, float[][][] pa, float[][][] pb) 
  {
    int n1 = pa[0][0].length;
    int n2 = pa[0].length;
    int nc = pa.length;
    float[][][] pc = new float[nc][n2][n1];
    for (int ic=0; ic<nc; ++ic) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          pc[ic][i2][i1] = (i1<depth)?pa[ic][i2][i1]:pb[ic][i2][i1];
        }
      }
    }
    return pc;
  }

  private static float[][][][] combine(
    float depth, float[][][][] pa, float[][][][] pb) 
  {
    int n1 = pa[0][0][0].length;
    int n2 = pa[0][0].length;
    int n3 = pa[0].length;
    int nc = pa.length;
    float[][][][] pc = new float[nc][n3][n2][n1];
    for (int ic=0; ic<nc; ++ic) {
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          pc[ic][i3][i2][i1] = (i1<depth)?pa[ic][i3][i2][i1]:pb[ic][i3][i2][i1];
        }
      }}
    }
    return pc;
  }


  private static float[][][] addWavelet(double fpeak, float[][][] p) {
    double sigma = max(1.0,1.0/(2.0*PI*fpeak));
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][] p0 = p[0];
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    float[][] q = copy(p0);
    float[][] q1 = new float[n2][n1];
    float[][] q2 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    for (int id=0; id<2; ++id) { // 2nd directional derivative of Gaussian
      rgf.apply10(q,q1);
      rgf.apply01(q,q2);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          q[i2][i1] = p1[i2][i1]*q1[i2][i1]+p2[i2][i1]*q2[i2][i1];
        }
      }
    }
    q = mul(q,-1.0f/rms(q)); // negate for Ricker wavelet
    return new float[][][]{q,p1,p2};
  }

  private static float[][][][] addWavelet(double fpeak, float[][][][] p) {
    double sigma = max(1.0,1.0/(2.0*PI*fpeak));
    int n1 = p[0][0][0].length;
    int n2 = p[0][0].length;
    int n3 = p[0].length;
    float[][][] p0 = p[0];
    float[][][] p1 = p[1];
    float[][][] p2 = p[2];
    float[][][] p3 = p[3];
    float[][][] q = copy(p0);
    float[][][] q1 = new float[n3][n2][n1];
    float[][][] q2 = new float[n3][n2][n1];
    float[][][] q3 = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    for (int id=0; id<2; ++id) { // 2nd directional derivative of Gaussian
      rgf.apply100(q,q1);
      rgf.apply010(q,q2);
      rgf.apply001(q,q3);
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            q[i3][i2][i1] = p1[i3][i2][i1]*q1[i3][i2][i1] +
                            p2[i3][i2][i1]*q2[i3][i2][i1] +
                            p3[i3][i2][i1]*q3[i3][i2][i1];
          }
        }
      }
    }
    q = mul(q,-1.0f/rms(q)); // negate for Ricker wavelet
    return new float[][][][]{q,p1,p2,p3};
  }

  private static float rms(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    double sum = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fi = f[i2][i1];
        sum += fi*fi;
      }
    }
    return (float)sqrt(sum/n1/n2);
  }
  private static float rms(float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    double sum = 0.0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fi = f[i3][i2][i1];
          sum += fi*fi;
        }
      }
    }
    return (float)sqrt(sum/n1/n2/n3);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
