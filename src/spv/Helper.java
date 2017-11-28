/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package spv;

import util.*;
import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Generates 2D fake training data for DeepFault.
 * <em>
 * Jacobians of functions used in folding and faulting have been implemented
 * but not tested. Therefore, beware of errors in calculated slopes p2 and p3.
 * </em>
 * @author Xinming Wu, University of Texas at Austin.
 * @version 2017.08.01
 */
public class Helper {

  public float[][][] convert(byte[][][] ep) {
    int n3 = ep.length;
    int n2 = ep[0].length;
    int n1 = ep[0][0].length;
    float[][][] ec = new float[n1][n2][n3];
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      ec[i1][i2][i3] = (float)ep[i3][i2][i1];
    }}}
    return ec;
  }

  public float[][][] fill(float[][][] ep) {
    int n3 = ep.length;
    int n2 = ep[0].length;
    int n1 = ep[0][0].length;
    float[][][] ef = new float[n3][n2][n1];
    for (int i1=0; i1<n1; i1++) {
    for (int i3=1; i3<n3-1; i3++) {
    for (int i2=1; i2<n2-1; i2++) {
      float epa = ep[i3-1][i2][i1];
      float epb = ep[i3+1][i2][i1];
      float epl = ep[i3][i2-1][i1];
      float epr = ep[i3][i2+1][i1];
      if(epa>0f&&epb>0f) ef[i3][i2][i1] = 1f;
      if(epl>0f&&epr>0f) ef[i3][i2][i1] = 1f;
    }}}
    return ef;
  }


}
