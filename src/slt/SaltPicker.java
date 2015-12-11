package slt;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

import ad.*;

/**
 * Extract salt boundaries from a salt likelihood image. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.10
 */

public class SaltPicker {

  public float[][] polarCoordinates(float c1, float c2, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float d1 = c1*c1+c2*c2;
    float d2 = (c1-n1+1)*(c1-n1+1)+c2*c2;
    float d3 = c1*c1+(c2-n2+1)*(c2-n2+1);
    float d4 = (c1-n1+1)*(c1-n1+1)+(c2-n2+1)*(c2-n2+1);
    float dm = max(d1,d2);
    dm = max(dm,d3);
    dm = max(dm,d4);
    int nr = round(sqrt(dm));
    float[][] fr = new float[360][nr*2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float xc1 = i1-c1;
      float xc2 = i2-c2;
      float xcs = sqrt(xc1*xc1+xc2*xc2);
      int ri = round(xcs);if(ri>=nr){ri=nr-1;}
    }}
    return fr;
  }
}
