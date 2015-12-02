package aii;

import edu.mines.jtk.dsp.FftComplex;


public class PredictionErrorFilter3 {
  private final float[] _cR;
  private final float[] _cI;
  
  private float[] _fR;
  private float[] _fI;

  // No test for cR and cI being of equal length
  public PredictionErrorFilter3( float[] cR, float[] cI ) {
    _cR = cR;
    _cI = cI;
  }

  public final float[] getReal() {
      return _fR;
  }
  public final float[] getImag() {
      return _fI;
  }

  public final void compute( int nout ) {
    
      // Set up filter elements (real and imaginary)
    _fR = new float[nout];
    _fI = new float[nout];
    
    // Set up power spectrum (length power of 2) and initialize to 1
    int np2 = (int)(Math.log((_cR.length))/Math.log(2.0));
    np2 = (int) Math.pow(2.0, (double)np2);
    if ( np2 < _cR.length ) np2 *= 2;
    np2 = FftComplex.nfftFast( np2 );
    double[] ps = new double[np2];
    for (int i = 0; i < ps.length; i++) {
      ps[i] = 1.0;
    }

    // Set up work arrays for input (odd is a copy, even is shifted by one)   
    float[] oR = new float[_cR.length];  // Mathematical indices 1, 3, 5, ...
    float[] oI = new float[_cR.length];
    float[] eR = new float[_cR.length];  // Mathematical indices 2, 4, 6, ...
    float[] eI = new float[_cR.length];
    int ii;
    for (int i = 0; i < _cR.length-1; i++) {   
        ii = i + 1;
      oR[i] = _cR[i];   oI[i] = _cI[i];
      eR[i] = _cR[ii];  eI[i] = _cI[ii];
    }   

    // Set up work arrays
    double p = 1.0;   
    float[] pf  = new float[2*np2];
    float[] wkR = new float[nout];
    float[] wkI = new float[nout];

    //Loop through specified number of samples in output spectrum
    FftComplex cfft = new FftComplex( np2 );
    for (int j = 0; j < nout; j++) {
      // Forward prediction of current element f[j] using odd and even elements
      // predicted to this point
    
      // Compute summed squared products & cross products of odd & even elements
      // psum = summation( (Real(i+1)+Real(i))**2 + (Imag(i+1)+Imag(i))**2 )
      // pdif = summation( (Real(i+1)-Real(i))**2 + (Imag(i+1)-Imag(i))**2 )
      // xsum = summation( (Real(i+1)+Imag(i))**2 + (Imag(i+1)-Real(i))**2 )
      // xdif = summation( (Real(i+1)-Imag(i))**2 + (Imag(i+1)+Real(i))**2 )
      double psum = 0.0, pdif = 0.0, xsum = 0.0, xdif = 0.0, dr, di;
      for (int i = 0; i < _cR.length-j-1; i++) {     
        dr = eR[i] + oR[i];  di = eI[i] + oI[i];  psum += (dr*dr + di*di);
        dr = eR[i] - oR[i];  di = eI[i] - oI[i];  pdif += (dr*dr + di*di);
        dr = eR[i] + oI[i];  di = eI[i] - oR[i];  xsum += (dr*dr + di*di);
        dr = eR[i] - oI[i];  di = eI[i] + oR[i];  xdif += (dr*dr + di*di);
      }

      // Compute element - Real = Squared sum of difference products minus
      //                          Squared sum of summed products
      //                   Imag = Squared sum of summed cross-products minus
      //                          Squared sum of difference cross-products
      // Normalize with sum of Squared sum of summed products and 
      //                       Squared sum of difference products     
      double d = psum + pdif;
      if ( d < 1.0e-10 ) {
        _fR[j] = _fI[j] = 0.0f;
      }
      else {
        _fR[j] = (float)((pdif-psum)/d);
        _fI[j] = (float)((xsum-xdif)/d);
//        _fI[j] = (float)((xdif-xsum)/d);
      }
      
      // Accumulate power spectrum normalizer product with (1.0-norm of f[j])
      p *= (1.0f - (_fR[j]*_fR[j]+_fI[j]*_fI[j]));

      // Update previous elements with new element, work[j]
      // Add to each f[i] the mult. of f[j] with the conjugate of f[j-i]
      for (int i = 0; i < j; i++) {
        wkR[i] = _fR[i];
        wkI[i] = _fI[i];
      }
      int jm1 = j - 1;
      for (int i = 0; i < j; i++) {
          ii = jm1 - i;
        _fR[i] += ( wkR[ii]*_fR[j] + wkI[ii]*_fI[j]);
        _fI[i] += (-wkI[ii]*_fR[j] + wkR[ii]*_fI[j]);
      }

      // Forward transform current solution
      
      pf[0] = 1.0f;  pf[1] = 0.0f;
      int jj = 2; 
      for (int i = 1; i <= j+1; i++) {
        ii = i - 1;
        pf[jj++] = _fR[ii];  
        pf[jj++] = _fI[ii];
      }
      for (int i = jj; i < 2*np2; i++) {
        pf[i] = 0.0f;
      }
      
      cfft.complexToComplex( 1, pf, pf );
      float[] pfR = new float[np2];
      float[] pfI = new float[np2];
      ii = 0;
      for (int i = 0; i < np2; i++) {
        pfR[i] = pf[ii++];
        pfI[i] = pf[ii++];
      }
      
      // Update power spectrum with normalized squared product of absolute 
      // values of transform
      float a;
      for (int i = 0; i < np2; i++) {
        if ( pfR[i] == 0.0f ) 
            a = Math.abs( pfI[i] );
        else if ( pfI[i] == 0.0f ) 
            a = Math.abs( pfR[i] );
        else {
        float aR = Math.abs( pfR[i] );
        float aI = Math.abs( pfI[i] );
        if ( aR > aI ) {
          float t = aI / aR;
          a = aR * (float)Math.sqrt( 1.0f + t*t );
        }
        else {
          float t = aR / aI;
          a = aI * (float)Math.sqrt( 1.0f + t*t );
        }
        }
      ps[i] += ( a*a / p );
    }

      // Update odd and even work arrays
      //   Odd  = odd + even * conjugate of f[j]
      //   Even = even[i+1] + odd[i+1] * f[j]
      if ( j < nout-1 ) {
        for (int i = 0; i < _cR.length-j-2; i++) {
            ii = i + 1;
          oR[i] += ( _fR[j]*eR[i] + _fI[j]*eI[i]);
          oI[i] += (-_fI[j]*eR[i] + _fR[j]*eI[i]); 
          eR[i] = eR[ii] + (_fR[j]*oR[ii] - _fI[j]*oI[ii]);
          eI[i] = eI[ii] + (_fI[j]*oR[ii] + _fR[j]*oI[ii]);
        }
      }
    }

    // Invert power spectrum and inverse fourier transform (complex to complex)
    ii = 0;
    for (int i = 0; i < np2; i++) {
      pf[ii++] = (float)(1.0/ps[i]);  
      pf[ii++] = 0.0f;
    }
    cfft.complexToComplex( -1, pf, pf );
    cfft.scale( np2, pf );   
    
    // Stabilization is similar to adding small noise to diagonal in time domain
    float[] pfR = new float[np2];
    float[] pfI = new float[np2];
    ii = 0;
    for (int i = 0; i < np2; i++) {
      pfR[i] = pf[ii++];
      pfI[i] = pf[ii++];
    }
    p = 1.000001 * pfR[0];
    for (int j = 0; j < nout; j++) {
        _fR[j] = _fI[j] = 0.0f;
      for (int i = 0; i < j; i++) {
        int jmi = j - i;
        _fR[j] += (_fR[i]*pfR[jmi] - _fI[i]*pfI[jmi]);
        _fI[j] += (_fI[i]*pfR[jmi] + _fR[i]*pfI[jmi]); 
      }
      _fR[j] = (_fR[j] + pfR[j+1]) / (float)(-p);
      _fI[j] = (_fI[j] + pfI[j+1]) / (float)(-p);

      for (int i = 0; i < j; i++) {
        wkR[i] =  _fR[i];  wkI[i] = -_fI[i];
      }
      int jm1 = j - 1;
      for (int i = 0; i < j; i++) {
          ii = jm1 - i;
        _fR[i] += (wkR[ii]*_fR[j] - wkI[ii]*_fI[j]);
        _fI[i] += (wkI[ii]*_fR[j] + wkR[ii]*_fI[j]);
      }
      
      float a;
      if ( _fR[j] == 0.0f ) 
        a = Math.abs( _fI[j] );
      else if ( _fI[j] == 0.0f ) 
        a = Math.abs( _fR[j] );
      else {
        float aR = Math.abs( _fR[j] );
        float aI = Math.abs( _fI[j] );
        if ( aR > aI ) {
          float t = aI / aR;
          a = aR * (float)Math.sqrt( 1.0f + t*t );
        }
        else {
          float t = aR / aI;
          a = aI * (float)Math.sqrt( 1.0f + t*t );
        }
      }
      p *= (1.0 - a*a);
    }
  }
}
