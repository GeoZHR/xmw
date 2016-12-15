/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package acm;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.mesh.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Active contour
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */

public class ActiveContour3 {


  public SnakeMesh initialSnake(float c1, float c2, float c3, float r) {
    float dt = 5f;
    float dp = 5f;
    float ft = 0.f;
    float fp = 0.f;
    int nt = round((180f-ft)/dt);
    int np = round((360f-fp)/dt);
    Sampling st = new Sampling(nt,dt,ft);
    Sampling sp = new Sampling(np,dp,fp);
    SnakeMesh sm = new SnakeMesh();
    SnakeMesh.Node[] nds = new SnakeMesh.Node[6];
    nds[0] = new SnakeMesh.Node(c3-r,c2,c1);
    nds[1] = new SnakeMesh.Node(c3+r,c2,c1);
    nds[2] = new SnakeMesh.Node(c3,c2-r,c1);
    nds[3] = new SnakeMesh.Node(c3,c2+r,c1);
    nds[4] = new SnakeMesh.Node(c3,c2,c1-r);
    nds[5] = new SnakeMesh.Node(c3,c2,c1+r);
    sm.addNodes(nds);
    for (int it=0; it<nt; ++it) {
      float ti = (float)toRadians(st.getValue(it));
    for (int ip=0; ip<np; ++ip) {
      float pi = (float)toRadians(sp.getValue(ip));
      float x = c3+r*sin(ti)*cos(pi);
      float y = c2+r*sin(ti)*sin(pi);
      float z = c1+r*cos(ti);
      SnakeMesh.Node nd = new SnakeMesh.Node(x,y,z);
      sm.addNode(nd);
    }}
    /*
    SnakeMesh.Node nr = sm.findNodeNearest(c3,c2,c1+r);
    sm.removeNode(nr);
    SnakeMesh.Node nd = new SnakeMesh.Node(c3,c2,c1+r+20);
    sm.addNode(nd);
    */
    return sm;
  }

  public void fakeField(float c1, float c2, float c3, float r, float[][][] f) {
    int b1 = round(c1-r);
    int e1 = round(c1+r);
    int b2 = round(c2-r);
    int e2 = round(c2+r);
    int b3 = round(c3-r);
    int e3 = round(c3+r);
    for (int i3=b3;i3<=e3;i3++){
    for (int i2=b2;i2<=e2;i2++){
      f[i3][i2][b1] = 1f;
      f[i3][i2][e1] = 1f;
    }}
    for (int i2=b2;i2<=e2;i2++){
    for (int i1=b1;i1<=e1;i1++){
      f[b3][i2][i1] = 1f;
      f[e3][i2][i1] = 1f;
    }}
    for (int i3=b3;i3<=e3;i3++){
    for (int i1=b1;i1<=e1;i1++){
      f[i3][b2][i1] = 1f;
      f[i3][e2][i1] = 1f;
    }}
  }



  public void updateSnake(
    int niter, SnakeMesh sm, float[][][] u1, float[][][] u2, float[][][] u3) {
    for (int it=0; it<niter; it++) {
      updateSnake(sm,u1,u2,u3);
    }
  }

  public void updateSnake(
    SnakeMesh sm, float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    int nd = sm.countNodes();
    SnakeMesh.Node[] ns = new SnakeMesh.Node[nd];
    SnakeMesh.NodeIterator ni = sm.getNodes();
    int k = 0;
    while(ni.hasNext()) {ns[k++] = ni.next();}
    float[] xm =  new float[3];
    float[] nv =  new float[3];
    ArrayList<SnakeMesh.Node> nl = new ArrayList<SnakeMesh.Node>();
    for (int id=0; id<nd; ++id) {
      SnakeMesh.Node ndi = ns[id];
      float x = ndi.x();
      float y = ndi.y();
      float z = ndi.z();
      int ix = round(x);
      int iy = round(y);
      int iz = round(z);
      iz = max(iz,0);
      iy = max(iy,0);
      iz = max(ix,0);
      ix = min(ix,n3-1);
      iy = min(iy,n2-1);
      iz = min(iz,n1-1);
      ndi.meanPosition(xm);
      ndi.normalVector(nv);
      float mx = xm[0];
      float my = xm[1];
      float mz = xm[2];
      float vx = nv[0];
      float vy = nv[1];
      float vz = nv[2];
      float fmx = mx-x;
      float fmy = my-y;
      float fmz = mz-z;
      float fms = fmx*vx+fmy*vy+fmz*vz;
      vx *= fms;
      vy *= fms;
      vz *= fms;
      float ftx = fmx-vx;
      float fty = fmy-vy;
      float ftz = fmz-vz;
      float uxi = u3[ix][iy][iz];
      float uyi = u2[ix][iy][iz];
      float uzi = u1[ix][iy][iz];
      float uvs = uxi*nv[0]+uyi*nv[1]+uzi*nv[2];
      //float xui = x+0.1f*vx+0.1f*ftx+nv[0]*uvs*10f;
      //float yui = y+0.1f*vy+0.1f*fty+nv[1]*uvs*10f;
      //float zui = z+0.1f*vz+0.1f*ftz+nv[2]*uvs*10f;
      float xui = x+ftx;//+uxi;
      float yui = y+fty;//+uyi;
      float zui = z+ftz;//+uzi;
      if(Float.isNaN(xui)||Float.isNaN(yui)||Float.isNaN(zui)){continue;}
      SnakeMesh.Node nu = new SnakeMesh.Node(xui,yui,zui);
      nl.add(nu);
    }
    int np = nl.size();
    SnakeMesh.Node[] nn = new SnakeMesh.Node[np];
    for (int ip=0; ip<np; ++ip)
      nn[ip] = nl.get(ip);
    sm.removeNodes(ns);
    sm.addNodes(nn);
  }
}
