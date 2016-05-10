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
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.14
 */

public class RosePlot {

  public void roseInDepths (float[][] fps, int nbin) {
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
    System.out.println("test1");
    float[] bins = applyHistogram(nbin,phi);
    System.out.println("test2");
    applyForRose(bins);
  }

  private void applyForRose(int nrow, int ncol,float[] bins) {
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
    pf.setFontSizeForPrint(8,240);
    pf.setVisible(true);
    pf.setSize(750,800);
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
    pf.setFontSizeForPrint(8,240);
    pf.setVisible(true);
    pf.setSize(750,800);
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
      float aci = (float)(0.5*phi/DBL_PI);
      Color rgb = Color.getHSBColor(aci,1f,1f);//getNextColor();
      float[][] rp1 = makeRadialPoints(rdi,phi);
      float[][] rp2 = makeRadialPoints(rdi,phi+dp);
      float dx = rp2[0][1] - rp1[0][1];
      float dy = rp2[1][1] - rp1[1][1];
      float ds = sqrt(dx*dx+dy*dy);
      float dd = 0.001f;
      float dk = dd*4f;
      float de = ds-dk;
      dx /= ds;
      dy /= ds;
      float di = dk;
      for ( ; di<=de; di+=dd) {
        float xi = rp1[0][1]+dx*di;
        float yi = rp1[1][1]+dy*di;
        float[] xs = new float[]{0f,xi};
        float[] ys = new float[]{0f,yi};
        PointsView pvr = pp.addPoints(xs,ys);
        if ((di-dk)<=dd||(di+dd)>=de) {
          pvr.setLineColor(Color.BLACK);
          pvr.setLineWidth(4f);
        } else { 
          pvr.setLineColor(rgb);
          pvr.setLineWidth(4f);
        }
      }
      float[] rx3 = new float[]{rp1[0][1]+(dk+dd)*dx,rp1[0][1]+(di-dd)*dx};
      float[] ry3 = new float[]{rp1[1][1]+(dk+dd)*dy,rp1[1][1]+(di-dd)*dy};
      float[][] rp3 = new float[][]{rx3,ry3};
      PointsView pvr3 = pp.addPoints(rp3[0],rp3[1]);
      pvr3.setLineColor(Color.BLACK);
      pvr3.setLineWidth(4.f);

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
