# Split Step Fourier Method 

Simulate the propagation of optical pulses in optical fiber with split step Fourier method (SSFM).
The propagation is governed by the nonlinear Schrödinger equation (NLSE).
We cannot solve the NLSE analytically and, thus, have to rely on numerical methods as in this repo.

### Things to do now
- Simulations with single polarization in multiple-span link
* How to simulate EDFA (add noise)?
* How to simulate CDF? Nonlinearity, dispersion, loss, etc.
* What DSP should be included? Currently chromatic dispersion (D) is compensated, but third order dispersion (S) is not. Should I consider other impairments?
  - $$a=b$$
* For the question above, I guess nonlinearity should not be compensated.
