# Split Step Fourier Method 

Simulate the propagation of optical pulses in optical fiber with split step Fourier method (SSFM).
The propagation is governed by the nonlinear Schrödinger equation (NLSE).
We cannot solve the NLSE analytically and, thus, have to rely on numerical methods as in this repo.

### To-do
- Simulations with single polarization in multiple-span link
* How to simulate EDFA (add noise)?
* How to simulate CDF? Nonlinearity, dispersion, loss, etc.
* ~~Currently chromatic dispersion (D) is compensated, but third order dispersion (S) is not. Should I consider other impairments?~~
  - [The lecture notes from the Fiber Optical Communication course in Chalmers](papers/dispersion-lecture-notes.pdf) introduce the dispersion in fiber. 
  - [A document on chromatic dispersion in optical fibers](papers/dispersion-general.pdf) summarizes how to computer third order dispersion coefficient.
* What DSP should be included? I guess nonlinearity should not be compensated.
  - No, only compensate for chromatic dispersion. The resultant constellation diagram of 16QAM is a bit rotated due to nonlinearity.
  - Back-propagation to the received signals can generate perfect constellation diagrams for all the WDM channels (including legacy OOK and 16QAM), so the SSF simulation should be correct.
  - In single polarization case, dispersion is enough for my simulations.


### Questions
- Spectral channel shape is RRC with 0.2 roll-off. The filter is normal RRC or square-root RRC?
- No shaping means Gaussian pulse in the time domain? Or rectangular shape in the spectrum domain?
- DCF: ideal FBG without loss, nonlienarity, etc.?
- Target pre-FEC Q=7.3 dB seems not very useful in simulations. It just gives some kind of threshold. But we are interested in the nonlinear noise suffered by the CUT in this transmission scenario, which should be represented by a numerical model or formula.
- If amplifiers are ideal noise-less, how is it possible to add equivalent noise at the end? The end of the whole transmission, or the end of each SSMF span? What is the noise figure, e.g., 5 dB?
