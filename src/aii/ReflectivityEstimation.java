package aii;

import edu.mines.jtk.dsp.FftReal;


public class ReflectivityEstimation {
  private float[] _seis;
  private float   _fextra, _spike, _dt;

  public ReflectivityEstimation() {
  }

  public ReflectivityEstimation( float[] seis, float fhigh, 
                                  float spike,  float dt ) {
    setParms( seis, fhigh, spike, dt);
  }
  
  public ReflectivityEstimation( float[] seis, float fhigh, 
                                 float spike,  float dt, int ns ) {
    setParms( seis, fhigh, spike, dt);
  }

  public ReflectivityEstimation( float fhigh, float spike, float dt) {
    setParms( fhigh, spike, dt);
  }
  
  public final void setParms( float fhigh, float spike, float dt) {
    _fextra = fhigh;
    _spike  = spike;
    _dt = dt;
  }
  
  public final void setParms( float[] seis, float fextra, 
                              float spike,  float dt) {
    setParms( fextra, spike, dt);
    setTrace( seis );
  }
  
  public void setTrace( float[] seis ) {
    _seis = seis;
  }

  public float[] compute() {

    // Determine non-zero part of trace; return if all zeros
    int itop = -1, ibot = -1;
    float smax = 0.0f;
    for (int i = 0; i < _seis.length; i++) {
      if ( _seis[i] != 0.0f ) {
        smax = Math.max( Math.abs(_seis[i]), smax);
        if ( itop < 0 ) 
          itop = i;
        else 
          ibot = i;
      }
    }
    if ( itop == -1 || ibot == -1 ) return _seis;
    if ( ibot < 0 ) return _seis;    
    smax *= 5.0f;
    
    int nt = FftReal.nfftFast( ibot-itop+1 );
    float[] values = new float[nt];
    for (int i = itop, ii = 0; i <= ibot; i++, ii++) {
      values[ii] = _seis[i] / smax;
    }
    
    // Generate complex transform of input trace
    FftReal fft = new FftReal( nt );
    float[] cvalues = new float[nt+2];
    fft.realToComplex( 1, values, cvalues );
    int nf = nt/2 + 1;
    float[] cR = new float[nf];
    float[] cI = new float[nf];
    for (int i = 0, ii = 0; i < cvalues.length-1; i+=2, ii++) {
      cR[ii] = cvalues[i];
      cI[ii] = cvalues[i+1];
    }

    // Truncate complex array of input through ih and 
    //   pad with conjugate of input ih to 2ih
    // ih - frequency index corresponding to extrapolation frequency
//    int ih = (int) (2.0*_dt*_fextra*(nf-1) + 1);
    int ih = (int) (2.0*_dt*_fextra*(nf-1) );
    float[] ccR = new float[2*ih+1];
    float[] ccI = new float[2*ih+1];
    for (int i = 1, i1 = ih+1, i2 = ih-1; i <= ih; i++, i1++, i2--) {
      ccR[i1] = cR[i];  ccI[i1] =  cI[i];
      ccR[i2] = cR[i];  ccI[i2] = -cI[i];
    }
    ccR[ih] = cR[0];  ccI[ih] = cI[0];
    
    // Initialization Step 1 as defined by Walker and Ulrych (1983)
    // Generate prediction error filter corresponding to specified extrapolation
    // frequency
    PredictionErrorFilter3 pef = new PredictionErrorFilter3( ccR, ccI );
    pef.compute( (Math.min( (int)(0.7*(ih+1)), 50)) );
    float[] pefR = pef.getReal();
    float[] pefI = pef.getImag();

    // Apply pef to complex input array above extrapolation frequency
    for (int j = ih; j < nf; j++) {  // Changed from j = ih+1
      cR[j] = 0.0f;  cI[j] = 0.0f;
      int jm1 = j - 1;
      for (int i = 0; i < pefR.length; i++) {
        int  ii = jm1 - i ;
        cR[j] -= (cR[ii]*pefR[i] - cI[ii]*pefI[i]);
        cI[j] -= (cI[ii]*pefR[i] + cR[ii]*pefI[i]);
      }
    }
    cI[nf-1] = 0.0f;
    for (int i = 0; i < nf; i++) {
      cR[i] *= 2.0f;  cI[i] *= 2.0f;
    }

    // Copy complex array and pad to twice the length of input complex array 
    // Back transform to get rcs trace with pef applied
    int nt2 = FftReal.nfftFast( 2*nt );
    float[] s1 = new float[nt2];
    cvalues = new float[nt2+2];
    for (int i = 0, ii = 0; i < nf; i++, ii+=2) {
      cvalues[ii]   = cR[i];
      cvalues[ii+1] = cI[i];
    }

    fft = new FftReal( nt2 );
    fft.complexToReal( -1, cvalues, s1 );
    fft.scale( nt2, s1 );

    // Generate auto gain window values
    // Idea of gain application and subsequent removal to improve results 
    // credited to Pusey and Wilkinson of Chevron Oil Field Research Co.
    int ng = 201;  // This should be examined
    float[] g = new float[s1.length];
    int ngm1 = ng - 1;
    for (int j = 0; j < s1.length; j++) {
      float gmax = Math.abs( s1[j] );
      for (int i = 0; i < Math.min( s1.length-j, ngm1 ); i++) {
        gmax = Math.max( (ng-i)*Math.abs(s1[j+i])/ng, gmax);
      }
      for (int i = 0; i < Math.min( j-1, ngm1 ); i++) {
        gmax = Math.max( (ng-i)*Math.abs(s1[j-i])/ng, gmax );
      }
      g[j] = 1.0f / gmax;
    }
    
    // Apply auto gain to input
    float s1max = Math.abs(s1[0]);
    for (int i = 1; i < s1.length; i++) {
      s1[i] *= g[i];
      s1max = Math.max( Math.abs(s1[i]), s1max );
    }
    
    // Set up work arrays
    float[] s2 = new float[s1.length];
    float[] s3 = new float[s1.length];
    System.arraycopy( s1, 0, s2, 0, s1.length );
    System.arraycopy( s1, 0, s3, 0, s1.length );

    // Compute initial complex array for gained input
    fft.realToComplex( 1, s2, cvalues );
    float[] c1R = new float[cvalues.length/2];
    float[] c1I = new float[cvalues.length/2];
    for (int i = 0, ii = 0; i < cvalues.length-1; i+=2, ii++) {
      c1R[ii] = cvalues[i];
      c1I[ii] = cvalues[i+1];
    }

    // Set up iteration parameters   
    int iter  = 0;
    int itest = -100;
    int n1 = s2.length * (ih-1)/s1.length + 1;
    
//    float aa = 10.0f / s1max / _spike;  // Test numerator (e.g. 10)
    float aa = 5.0f / s1max / _spike;  // Test numerator (e.g. 10)
    aa = 1 / s1max;
    
    // Iterate extrapolation until minimum tolerance [change in L2/(L1*L1)]
    // Referred to as the EFMED method by Walker and Ulrych (1983)
    double tolerance = 5.0e-03;  // Test for convergence setting
    float  normLast  = 1.0e30f;
    boolean nextIteration = true;
    while ( nextIteration ) {
      iter++;
      
      // Compute L1 and L2 norms of the exponential transformation 
      //    (Ooe and Ulrych, 1979)
      // Transformation is a function of the degree of "spikiness" desired
      // Step 2
      float L1 = 0.0f;
      float L2 = 0.0f;
      for (int i = 0; i < s1.length; i++) {
        double b = aa * s1[i];
        s2[i] = 1.0f - (float) Math.exp( -b*b );
        L1 += s2[i];
        L2 += s2[i] * s2[i];
      }
      float norm = L2 / (L1*L1);
      
      // Apply exponential transformations (Equation 25) and reinitialize pad
      float a1 = L1 / L2;
      float a2 = 1.0f + a1;
      for (int i = 0; i < s1.length; i++) {
        s2[i] *= s3[i] * (a2 - a1*s2[i]);
      }

      // Transform updated rcs to frequency domain
      fft.realToComplex( 1, s2, cvalues );
      ccR = new float[cvalues.length/2];
      ccI = new float[cvalues.length/2];
      for (int i = 0, ii = 0; i < cvalues.length-1; i+=2, ii++) {
        ccR[ii] = cvalues[i];
        ccI[ii] = cvalues[i+1];
      }

      // Step 3 
      // Subtract modulated transformed array from initial within passband
      a1 = 1.0f / (L1+1.0f);
      a2 = 1.0f - a1;
      int ii = 0;
      for (int i = 0; i < n1; i++) {
        cvalues[ii++] = a2*c1R[i] + a1*ccR[i];
        cvalues[ii++] = a2*c1I[i] + a1*ccI[i];
      }
      for (int i = n1; i < cvalues.length/2; i++) {
        cvalues[ii++] = ccR[i];
        cvalues[ii++] = ccI[i];
      }

      // Compute iteration norm L2/(L1*L1), a modified varimax norm 
      //   (Ooe and Ulrych, 1979)
      // and compare with previous iteration
      if ( Math.abs(normLast-norm)/normLast < tolerance ) {
         if ( iter - itest == 1 ) nextIteration = false;
         if ( iter > 50000 ) nextIteration = false;
         itest = iter;
      }
      normLast = norm;

      // Back transform resulting reflectivity
      // Reinitialize pad
      fft.complexToReal( -1, cvalues, s3 );
      fft.scale( s3.length, s3 );
    }

    // Assign reflectivity to output array after gain removal
    // Pad back 0 to itop; ibot to original input length
    // Rescale back to input scaling
    float[] rcs = new float[_seis.length];
//    int ii = 0;
    for (int i = itop, ii = 0; i <= ibot; i++, ii+=2) {
      rcs[i] = s3[ii] / g[ii];
//      rcs[i] *= smax;
//      ii += 2;
    } 
    return rcs;
  }
  
  public float getFrequency( float[] trace, float dt, float amp ) {
    SpectralAnalysis spect = new SpectralAnalysis( trace, dt );
    spect.forward();
    float frequency = 0.0f;
    float df = 1.0f / (2.0f*dt) / (trace.length/2+1);
    float[] amplitude = spect.amplitude();
    for (int i = amplitude.length-1; i >= 0; i--) {
      if ( amplitude[i] >= amp ) {
        frequency = i*df;
        break;
      }
    }
    return frequency;
  }

/*
  public static void main( String[] args) {

    float fextra = 55.0f;
    float spike  =  1.0f;
    float dt     = 0.002f;

    TraceRegular seis = new TraceRegular();
    seis.readAsciiFile( "G:\\Data\\test\\trace_4430.dat", 501 );

    long start = System.currentTimeMillis();
    ReflectivityEstimation2 ref = new ReflectivityEstimation2( seis._values, 
                                 fextra, spike, 0.002f, 0 );
    float[] s = ref.compute();
    System.out.println("Float run time: "+(System.currentTimeMillis()-start));
//    for (int i = 0; i < s.length; i++) {
//      System.out.println( s[i] );
//    }
//    SpectralAnalysis2 spect = new SpectralAnalysis2( s, dt );
//    spect.forward();
//    float[] amp = spect.amplitude();
//    float   df  = spect.frequencyRate();
//    for (int i = 0; i < amp.length; i++) {
//      System.out.println((i*df)+"  "+amp[i]);
//    }
  }
*/
}
