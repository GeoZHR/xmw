package aii;

public class SpectralAnalysis {
  private float[]   _t, _amp = null, _faz = null;
  private float     _dt, _df, _avg;
  private Complex[] _f;
	
  public SpectralAnalysis() {}
	
  public SpectralAnalysis( float[] t, float dt ) {
    setSeries( t );
    _dt = dt;
  }

  public void setSeries( float[] t, float dt ) {
    setSeries( t );
    _dt = dt;
  }

  public void setSeries( float[] t ) {
    _t = t;
    _avg = 0.0f;
    for (int i = 0; i < _t.length; i++) {
      _avg += Math.abs( _t[i] );
    }
    _avg /= (float)_t.length;
  }
	
  public float[] series() {
    return _t;
  }
  public float[] amplitude() {
    return _amp;
  }
  public void setAmplitude( float[] amp ) {
    _amp = amp;
  }
  public float[] phase() {
    return _faz;
  }
  public void setPhase( float[] faz ) {
    _faz = faz;
  }
  public float timeRate() {
    return _dt;
  }
  public float frequencyRate() {
    return _df;
  }
  public void setFrequencyRate( float df ) {
    _df = df;
  }
	
  public void forward() {
    forward(true); 
  }
  
  public void forward(boolean normalize) {
	if ( _t == null ) return;
	RealFFT fft = new RealFFT();
	forward( fft , normalize );
  }
	
  public void forward( RealFFT fft ) {
    forward (fft, true);
  }
  
  public void forward( RealFFT fft, boolean normalize) {
    _f   = fft.forward( _t );
    _amp = new float[_f.length];
    _faz = new float[_f.length];
		
    if (_dt > 0.0f ) _df = 1.0f/(2.0f*_dt)/(float)(_f.length-1);

    float ampmax = 0.0f;
    float pih = (float)(Math.PI/2);
    for (int i = 0; i < _f.length; i++) {
      _amp[i] = (float)Math.sqrt( _f[i].norm() );
      ampmax = Math.max( _amp[i], ampmax );
      if ( _f[i]._r == 0.0f ) {
        if ( _f[i]._i > 0.0f )
          _faz[i] =  pih;
        else if (_f[i]._i < 0.0f )
          _faz[i] = -pih;
      }
      else
        _faz[i] = (float)Math.atan2( _f[i]._i, _f[i]._r );
    }

    if (normalize) {
      for (int i = 0; i < _amp.length; i++) {
        _amp[i] /= ampmax;
      }
    }
    
  }
  
  public void forwardDB() {
    forward();
    if ( _amp == null ) return;
    
    for (int i = 0; i < _amp.length; i++) {
      _amp[i] = 20.0f * (float)Math.log10( (double)_amp[i] );
    }
  }
	
  public void inverse() {
    if ( _amp == null || _faz == null ) return;
    RealFFT fft = new RealFFT();
    inverse( fft );
  }
	
  public void inverse( RealFFT fft ) {
		
    // Form complex array from amplitude and phase
    _f = new Complex[_amp.length];
    for (int i = 0; i < _amp.length; i++) {
      _f[i] = new Complex( _amp[i] * (float)Math.cos( _faz[i] ),
                           _amp[i] * (float)Math.sin( _faz[i] ) );
    }

    // Perform inverse FFT
    _t = fft.inverse( _f );
  }
  
  public void inverseDB() {
    if ( _amp == null ) return;
	    
    for (int i = 0; i < _amp.length; i++) {
      _amp[i] = (float)Math.pow( 10.0, (double)_amp[i]/20.0 );
    }
    inverse();
  }
	
  public void smooth( int n ) {
    if ( n <= 1 ) return;
    if ( _t == null ) return;
    if ( _amp == null ) forward();
		
    // Compute triangular weighted smoother (one-sided)
    float[] s = new float[n];
    s[0] = 1.0f;
    for (int i = 1; i < n; i++) {
      s[i] = (float)(n-i) / (float)n;
    }

    // Apply smoother (two-sided)
    float[] amp = new float[_amp.length];
    for (int i = 0; i < _amp.length; i++) {
      amp[i] = _amp[i] * s[0];
      int ii = i + 1;
      int jj = 1;
      while ( ii < _amp.length && jj < n ) {
        amp[i] += _amp[ii++] * s[jj++];
      }
      ii = i - 1;
      jj = 1;
      while ( ii >= 0 && jj < n ) {
        amp[i] += _amp[ii--] * s[jj++];
      }
    }
    float ampmax = amp[0];
    _amp[0] = amp[0];
    for (int i = 1; i < _amp.length; i++) {
      ampmax = Math.max( amp[i], ampmax );
      _amp[i] = amp[i];
    }  

    for (int i = 0; i < _amp.length; i++) {
      _amp[i] /= ampmax;
    }

	}
	
}
