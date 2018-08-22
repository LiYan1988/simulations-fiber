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
  - No, only compensate for chromatic dispersion. The resultant constellation diagram of 16QAM is a bit rotated due to __XPM__.
  - Back-propagation to the received signals can generate perfect constellation diagrams for all the WDM channels (including legacy OOK and 16QAM), so the SSF simulation should be correct.
  - In single polarization case, dispersion is enough for my simulations.


