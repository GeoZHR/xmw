/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package crf;

import java.awt.*;
import java.util.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

import mef.*;

/**
 * Make a rosette-strike plot for faults/fault cells.
 * @author Xinming Wu, University of Texas at Austin.
 * @version 2016.05.09
 */

public class RosePlot {

  public void rotate(float phi, float[] fps) {
    int np = fps.length;
    for (int ip=0; ip<np; ++ip) {
      float fpi = fps[ip];
      if(fpi>=0.0f) {
        fpi += phi;
        if (fpi>=360f) fpi-=360f;
        fps[ip] = fpi;
      }
    }
  }

  public void convert(float[] fps) {
    int np = fps.length;
    for (int ip=0; ip<np; ++ip) {
      float fpi = fps[ip];
      if(fpi>180.0f) {
        fpi -= 180f;
        fps[ip] = fpi;
      }
    }
  }


  public void applyForRosePlots(
    float b1, float e1, int c2, int c3, int n2, int n3, 
    int nbin, float[][] fp, float[][] ob) 
  {
    int np = fp[0].length;
    int d2 = round(n2/c2);
    int d3 = round(n3/c3);
    int wx = min(round(1450f/c2),400);
    int wy = wx+12;
    int j3 = -1;
    for (int i3=0; i3<c3; ++i3) {
    //for (int i3=c3-1; i3>=0; --i3) {
      j3++;
      int b3 = i3*d3;
      int e3 = (i3+1)*d3;
      if(i3==(c3-1)){e3 = max(e3,n3);}
      System.out.println("b3="+b3);
      System.out.println("e3="+e3);
    for (int i2=0; i2<c2; ++i2) {
      int b2 = i2*d2;
      int e2 = (i2+1)*d2;
      if(i2==(c2-1)){e2 = max(e2,n2);}
      System.out.println("b2="+b2);
      System.out.println("e2="+e2);
      ArrayList<Float> fpa = new ArrayList<Float>();
      for (int ip=0; ip<np; ++ip) {
        int k2 = (int)fp[1][ip];
        int k3 = (int)fp[2][ip];
        if(k2<b2 || k2>e2) {continue;}
        if(k3<b3 || k3>e3) {continue;}
        float obi = ob[k3][k2];
        int k1 = (int)fp[0][ip];
        if(k1>=(b1+obi) && k1<=(e1+obi)) {
        //if(k1>=b1 && k1<=e1) {
          float fpi = fp[3][ip];
          fpa.add(fpi);
        }
      }
      int npk = fpa.size();
      System.out.println("npk="+npk);
      float[] fpk = new float[npk];
      for (int ik=0; ik<npk; ++ik) {
        fpk[ik] = fpa.get(ik);
      }
      fpa.clear();
      int x = i2*wx;
      int y = j3*wy;
      if(j3>0) {y+=20;}
      rose(x, y, wx, wy,  fpk, nbin);
    }}
  }

  public void rose(float[][][] fp, int nbin) {
    int n3 = fp.length;
    int n2 = fp[0].length;
    int n1 = fp[0][0].length;
    ArrayList<Float> fpa = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fpi = fp[i3][i2][i1];
      if (fpi>=0.0f) {fpa.add(fpi);}
    }}}
    int ip = 0;
    int np = fpa.size();
    float[] fps = new float[np];
    for (float fpi:fpa) {
      fps[ip] = fpi;
      ip++;
    }
    rose(fps,nbin);
  }

  public float[][] faultPoints(float[][][] fp) {
    int n3 = fp.length;
    int n2 = fp[0].length;
    int n1 = fp[0][0].length;
    ArrayList<Float> fpa = new ArrayList<Float>();
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    ArrayList<Float> x3a = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fpi = fp[i3][i2][i1];
      if (fpi>=0.0f) {
        fpa.add(fpi);
        x1a.add((float)i1);
        x2a.add((float)i2);
        x3a.add((float)i3);
      }
    }}}
    int np = fpa.size();
    float[][] fps = new float[4][np];
    for (int ip=0; ip<np; ++ip) {
      fps[0][ip] = x1a.get(ip);
      fps[1][ip] = x2a.get(ip);
      fps[2][ip] = x3a.get(ip);
      fps[3][ip] = fpa.get(ip);
    }
    return fps;
  }

  public void rose(FaultSkin[] skins, int nbin) {
    FaultCell[] cells = FaultSkin.getCells(skins);
    int nc = cells.length;
    float[] fp = new float[nc];
    for (int ic=0; ic<nc; ++ic)
      fp[ic] = cells[ic].getFp();
    rose(fp,nbin);
  }

  public void rose(float[] phi, int nbin) {
    float[] bins = applyHistogram(nbin,phi);
    applyForRose(bins);
  }

  public void rose(int x, int y, int wx, int wy, float[] phi, int nbin) {
    float[] bins = applyHistogram(nbin,phi);
    applyForRose(x,y,wx,wy,bins);
  }



  private void applyForRose(int x, int y, int wx, int wy, float[] bins) {
    PlotPanel.AxesPlacement axes = PlotPanel.AxesPlacement.NONE;
    //PlotPanel.Orientation orient = PlotPanel.Orientation.X1DOWN_X2RIGHT;
    PlotPanel.Orientation orient = PlotPanel.Orientation.X1RIGHT_X2UP;
    PlotPanel pp = new PlotPanel(1,1,orient,axes);
    for (float rdi=0.2f; rdi<=1.0f; rdi+=0.2f) {
      float[][] cps = makeCirclePoints(1000,rdi);
      PointsView pvc = pp.addPoints(cps[0],cps[1]);
      pvc.setLineStyle(PointsView.Line.DASH);
      if(rdi==1.0f) {
        pvc.setLineColor(Color.RED);
        pvc.setLineWidth(4.f);
      } else {
        pvc.setLineColor(Color.BLACK);
        pvc.setLineWidth(3.f);
      }
    }
    addRadials(pp,Color.BLACK);
    //addBins(pp,bins);
    float rc = addBinsF(pp,bins);
    pp.setTitle("Red cycle: "+rc*100f+"%");
    PlotFrame pf = new PlotFrame(pp);
    pf.setFontSizeForPrint(10,480);
    pf.setVisible(true);
    pf.setSize(wx,wy);
    pf.setLocation(x,y);
  }


  private void applyForRose(float[] bins) {
    PlotPanel.AxesPlacement axes = PlotPanel.AxesPlacement.NONE;
    //PlotPanel.Orientation orient = PlotPanel.Orientation.X1DOWN_X2RIGHT;
    PlotPanel.Orientation orient = PlotPanel.Orientation.X1RIGHT_X2UP;
    PlotPanel pp = new PlotPanel(1,1,orient,axes);
    for (float rdi=0.2f; rdi<=1.0f; rdi+=0.2f) {
      float[][] cps = makeCirclePoints(1000,rdi);
      PointsView pvc = pp.addPoints(cps[0],cps[1]);
      pvc.setLineStyle(PointsView.Line.DASH);
      if(rdi==1.0f) {
        pvc.setLineColor(Color.RED);
        pvc.setLineWidth(4.f);
      } else {
        pvc.setLineColor(Color.BLACK);
        pvc.setLineWidth(3.f);
      }
    }
    addRadials(pp,Color.BLACK);
    //addBins(pp,bins);
    float rc = addBinsF(pp,bins);
    pp.setTitle("Red cycle: "+rc*100f+"%");
    PlotFrame pf = new PlotFrame(pp);
    pf.setFontSizeForPrint(10,480);
    pf.setVisible(true);
    pf.setSize(400,412);
  }

  private float[] applyHistogram(int nbins, float[] phi) {
    int np = phi.length;
    float ps = 1f/(float)np;
    float dp = 360f/(float)nbins;
    float[] bins = new float[nbins];
    for (int ip=0; ip<np; ++ip) {
      int ib = (int)(phi[ip]/dp);
      bins[ib] += ps;
    }
    return bins;
  }


  private void addRadials(PlotPanel pp, Color color) {
    int np = 8;
    float dp = (float)(2*DBL_PI/np);
    for (int ip=0; ip<np; ++ip) {
      float phi = dp*ip;
      float[][] rps = makeRadialPoints(1f,phi);
      PointsView pv = pp.addPoints(rps[0],rps[1]);
      pv.setLineStyle(PointsView.Line.DASH);
      pv.setLineColor(color);
      pv.setLineWidth(4.f);
    }
  }


  private void addBins(PlotPanel pp, float[] bins) {
    int nb = bins.length;
    float dp = (float)(2*DBL_PI/nb);
    float rmax = max(bins);
    mul(bins,0.8f/rmax,bins);
    for (int ib=0; ib<nb; ++ib) {
      float phi = ib*dp;
      float rdi = bins[ib];
      Color rgb = getNextColor();
      float[][] rp1 = makeRadialPoints(rdi,phi);
      float[][] rp2 = makeRadialPoints(rdi,phi+dp);
      /*
      float[][] cps = makeCirclePoints(1000,rdi);
      PointsView pvc = pp.addPoints(cps[0],cps[1]);
      pvc.setLineStyle(PointsView.Line.DASH);
      pvc.setLineColor(rgb);
      pvc.setLineWidth(3.f);
      */
      float[] rx3 = new float[]{rp1[0][1],rp2[0][1]};
      float[] ry3 = new float[]{rp1[1][1],rp2[1][1]};
      float[][] rp3 = new float[][]{rx3,ry3};
      PointsView pvr1 = pp.addPoints(rp1[0],rp1[1]);
      PointsView pvr2 = pp.addPoints(rp2[0],rp2[1]);
      PointsView pvr3 = pp.addPoints(rp3[0],rp3[1]);
      pvr1.setLineColor(rgb);
      pvr2.setLineColor(rgb);
      pvr3.setLineColor(rgb);
      pvr1.setLineWidth(4.f);
      pvr2.setLineWidth(4.f);
      pvr3.setLineWidth(4.f);
    }
  }

  private void addBinsB(PlotPanel pp, float[] bins) {
    int nb = bins.length;
    float dp = (float)(2*DBL_PI/nb);
    float rmax = max(bins);
    mul(bins,0.8f/rmax,bins);
    for (int ib=0; ib<nb; ++ib) {
      float phi = ib*dp;
      float rdi = bins[ib];
      Color rgb = getNextColor();
      float[][] rp1 = makeRadialPoints(rdi,phi);
      float[][] rp2 = makeRadialPoints(rdi,phi+dp);
      float dx = rp2[0][1] - rp1[0][1];
      float dy = rp2[1][1] - rp1[1][1];
      float ds = sqrt(dx*dx+dy*dy);
      for (float di=0f; di<=ds; di+=0.01f) {
        float xi = rp1[0][1]+dx*di/ds;
        float yi = rp1[1][1]+dy*di/ds;
        float[] xs = new float[]{0f,xi};
        float[] ys = new float[]{0f,yi};
        PointsView pvr = pp.addPoints(xs,ys);
        pvr.setLineColor(rgb);
        pvr.setLineWidth(4.f);
      }
      /*
      float[][] cps = makeCirclePoints(1000,rdi);
      PointsView pvc = pp.addPoints(cps[0],cps[1]);
      pvc.setLineStyle(PointsView.Line.DASH);
      pvc.setLineColor(rgb);
      pvc.setLineWidth(3.f);
      float[] rx3 = new float[]{rp1[0][1],rp2[0][1]};
      float[] ry3 = new float[]{rp1[1][1],rp2[1][1]};
      float[][] rp3 = new float[][]{rx3,ry3};
      PointsView pvr1 = pp.addPoints(rp1[0],rp1[1]);
      PointsView pvr2 = pp.addPoints(rp2[0],rp2[1]);
      PointsView pvr3 = pp.addPoints(rp3[0],rp3[1]);
      pvr1.setLineColor(rgb);
      pvr2.setLineColor(rgb);
      pvr3.setLineColor(rgb);
      pvr1.setLineWidth(4.f);
      pvr2.setLineWidth(4.f);
      pvr3.setLineWidth(4.f);
      */
    }
  }


  private float addBinsF(PlotPanel pp, float[] bins) {
    int nb = bins.length;
    float dp = (float)(2*DBL_PI/nb);
    float rmax = max(bins);
    float rc = round(rmax*120f)/100f;
    mul(bins,1.0f/rc,bins);
    for (int ib=0; ib<nb; ++ib) {
      float phi = ib*dp;
      float rdi = bins[ib];
      float aci = (float)(phi/DBL_PI);
      Color rgb = Color.getHSBColor(aci,1f,1f);
      float[][] rp1 = makeRadialPoints(rdi,phi);
      float[][] rp2 = makeRadialPoints(rdi,phi+dp);
      float dx = rp2[0][1] - rp1[0][1];
      float dy = rp2[1][1] - rp1[1][1];
      float ds = sqrt(dx*dx+dy*dy);
      float dd = 0.001f;
      float dk = dd*5f;
      float de = ds-dk;
      dx /= ds;
      dy /= ds;
      float xi,yi;
      float[] xs, ys;
      PointsView pvr;
      float di = dk;
      for ( ; di<de; di+=dd) {
        xi = rp1[0][1]+dx*di;
        yi = rp1[1][1]+dy*di;
        xs = new float[]{-xi,xi};
        ys = new float[]{-yi,yi};
        pvr = pp.addPoints(xs,ys);
        pvr.setLineColor(rgb);
        pvr.setLineWidth(3f);
      }
      // add black bounds
      xi = rp2[0][1];
      yi = rp2[1][1];
      xs = new float[]{-xi,xi};
      ys = new float[]{-yi,yi};
      pvr = pp.addPoints(xs,ys);
      pvr.setLineColor(Color.BLACK);
      pvr.setLineWidth(2f);

      xi = rp1[0][1];
      yi = rp1[1][1];
      xs = new float[]{-xi,xi};
      ys = new float[]{-yi,yi};
      pvr = pp.addPoints(xs,ys);
      pvr.setLineColor(Color.BLACK);
      pvr.setLineWidth(2f);

      float[][] rt1 = makeRadialPoints(rdi*1.01f,phi);
      float[][] rt2 = makeRadialPoints(rdi*1.01f,phi+dp);
      xs = new float[]{rt1[0][1],rt2[0][1]};
      ys = new float[]{rt1[1][1],rt2[1][1]};
      float[][] rp3 = new float[][]{xs,ys};
      pvr = pp.addPoints(rp3[0],rp3[1]);
      pvr.setLineColor(Color.BLACK);
      pvr.setLineWidth(2.f);
      pvr = pp.addPoints(mul(rp3[0],-1),mul(rp3[1],-1));
      pvr.setLineColor(Color.BLACK);
      pvr.setLineWidth(2.f);
    }
    return rc;
  }



  // Makes array of points on the unit circle in a plane.
  private float[][] makeCirclePoints(int nt, float r) {
    double dt = 2.0*DBL_PI/(nt-1);
    float[] x = new float[nt];
    float[] y = new float[nt];
    for (int it=0; it<nt; ++it) {
      float t = (float)(it*dt);
      x[it] = r*sin(t);
      y[it] = r*cos(t);
    }
    return new float[][]{x,y};
  }

  private float[][] makeRadialPoints(float r, float phi) {
    float[] x = new float[2];
    float[] y = new float[2];
    x[0] = 0.0f; x[1] = r*sin(phi); 
    y[0] = 0.0f; y[1] = r*cos(phi);
    return new float[][]{x,y};
  }

  private Color getNextColor() {
    float hue = (float)Math.random();
    float sat = (float)Math.random()*0.6f+0.4f;
    float bri = (float)Math.random()*0.3f+0.7f;
    int c = Color.HSBtoRGB(hue,sat,bri);
    int r = c&0x00ff0000;
    int g = c&0x0000ff00;
    int b = c&0x000000ff;
    return new Color(r|g|b);
  }
  
}
