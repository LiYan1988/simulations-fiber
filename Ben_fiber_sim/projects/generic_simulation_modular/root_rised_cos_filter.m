function u =root_rised_cos_filter(sig,fs,T,beta)


def _rrcos_pulseshaping_freq(sig, fs, T, beta):
    """
    Root-raised cosine filter in the spectral domain by multiplying the fft of the signal with the
    frequency response of the rrcos filter.

    Parameters
    ----------
    sig    : array_like
        input time distribution of the signal
    fs    : float
        sampling frequency of the signal
    T     : float
        width of the filter (typically this is the symbol period)
    beta  : float
        filter roll-off factor needs to be in range [0, 1]

    Returns
    -------
    sign_out : array_like
        filtered signal in time domain
    """
        
    f = np.fft.fftfreq(sig.shape[0])*fs
    nyq_fil = rrcos_freq(f, beta, T)
    nyq_fil /= nyq_fil.max()
    sig_f = np.fft.fft(sig)
    sig_out = np.fft.ifft(sig_f*nyq_fil)
    return sig_out