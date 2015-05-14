/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package util;


import edu.mines.jtk.sgl.*;

/**
 * make ellipsoids for displaying
 * @author Xinming Wu
 * @version 2015.04.13
 */
public class EllipsoidNode extends Node {
  public EllipsoidNode(
    float[] cx,float[] ux, float[] vx, float[] wx) {
    _cx = cx;
    _ux = ux;
    _vx = vx;
    _wx = wx;
  }
  public void draw(DrawContext dc) {
    float cx = _cx[2], cy = _cx[1], cz = _cx[0];
    float ux = _ux[2], uy = _ux[1], uz = _ux[0];
    float vx = _vx[2], vy = _vx[1], vz = _vx[0];
    float wx = _wx[2], wy = _wx[1], wz = _wx[0];
    _ellipsoid.draw(cx,cy,cz,ux,uy,uz,vx,vy,vz,wx,wy,wz);
  }
  public BoundingSphere computeBoundingSphere(boolean finite) {
    return _bs;
  }
  private float _thickness;
  private float[] _cx,_ux,_vx,_wx;
  private EllipsoidGlyph _ellipsoid = new EllipsoidGlyph();
  private BoundingSphere _bs = 
    new BoundingSphere(new BoundingBox(-1,-1,-1,1,1,1));
}
