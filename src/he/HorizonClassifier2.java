package he;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import vec.*;
import util.*;

/**
 * Extract a single seismic horizon curve by iterative 
 * seismic wavefome classification.
 * @author Xinming Wu
 * @version 2017.12.22
 */

public class HorizonClassifier2{

  public float[][] findPeaks() {
    return null;
  }
 
  private float semblance(float[] g, float[] f) {
    float sn = 0f;
    float sd = 0f;
    int n1 = g.length;
    for (int i1=0; i1<n1; ++i1) {

    }
    return sn/sd;
  }

  //A Shape-based Similarity Measure for Time Series Data 
  //with Ensemble Learning, 2013, proposed by Tetsuya Nakamura
  private float shapeBasedSimilarity(float[] g, float[] f) {
    float ss = 0;
    int n1 = g.length;
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i1];
      float fi = f[i1];
      float si = (1f+gi*fi)/sqrt((1f+gi*gi)*(1f+fi*fi)); 
      if(si>0) ss += si*_sv[i1];
    }
    return ss/n1;
  }

  private float[] _sh;
  private float[] _sv;
  private int _wh;
  private int _wv;
}
