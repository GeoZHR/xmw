/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package crf;

import java.awt.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

import mef.*;

/**
 * Make a rosette-strike plot for faults/fault cells.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.14
 */

public class RosePlot {

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


  private void applyForRose(float[] bins) {
    PlotPanel.AxesPlacement axes = PlotPanel.AxesPlacement.NONE;
    PlotPanel.Orientation orient = PlotPanel.Orientation.X1DOWN_X2RIGHT;
    PlotPanel pp = new PlotPanel(1,1,orient,axes);
    pp.setTitle("Rose-strike plot");
    float[][] cps = makeCirclePoints(1000,1);
    PointsView pvc = pp.addPoints(cps[0],cps[1]);
    pvc.setLineColor(Color.RED);
    pvc.setLineWidth(4.f);
    addRadials(pp,Color.BLACK);
    addBins(pp,bins);
    PlotFrame pf = new PlotFrame(pp);
    pf.setFontSizeForPrint(8,240);
    pf.setVisible(true);
    pf.setSize(750,800);
  }

  private float[] applyHistogram(int nbins, float[] phi) {
    int np = phi.length;
    float ps = 2f/(float)np;
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
    for (int ib=0; ib<nb; ++ib) {
      float phi = ib*dp;
      float rdi = bins[ib];
      Color rgb = getNextColor();
      float[][] rp1 = makeRadialPoints(rdi,phi);
      float[][] rp2 = makeRadialPoints(rdi,phi+dp);
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
    }
  }


  // Makes array of points on the unit circle in a plane.
  private float[][] makeCirclePoints(int nt, float r) {
    double dt = 2.0*DBL_PI/(nt-1);
    float[] x = new float[nt];
    float[] y = new float[nt];
    for (int it=0; it<nt; ++it) {
      float t = (float)(it*dt);
      x[it] = r*cos(t);
      y[it] = r*sin(t);
    }
    return new float[][]{x,y};
  }

  private float[][] makeRadialPoints(float r, float phi) {
    float[] x = new float[2];
    float[] y = new float[2];
    x[0] = 0.0f; x[1] = r*cos(phi); 
    y[0] = 0.0f; y[1] = r*sin(phi);
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
