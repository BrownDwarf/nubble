# Outline

1. What are the brightest known transiting-planet host stars? What will TESS find?

Table 3 of Rodriguez et al. 2017 shows the list of best 20 confirmed planets for transmission spectroscopy, for planets with radii < 5 R_earth.  We expect signal strengths in the ~100 parts per million range for targets in this list.

2. How small of a telescope can you measure such an exoplanet transit with?

- S/N versus telescope aperture for fixed J-mag
- S/N versus J-mag for fixed telescope aperture

3. How does this mission compare with existing or proposed missions?

- JWST  (slated for launch)
- CUTE (in development)
- FINESSE (proposed)
- Ground-based (Recent Giano paper, others)

4. What is the advantage of going to space?
- Stability
- Continuous monitoring (point and stare)
- Low/no background or telluric absorption

5. What is the nominal design of the spectrograph?
- Largest possible aperture that can fit in the cubesat
- Si Grism
- Custom MIT Lincoln Lab detector

6. What is the orbital configuration?
- Either Earth orbiting (TESS) or Earth trailing (Kepler, CUTE)
- Earth orbiting preferred for communications

7. What is the nominal operation of a cubesat with the proposed design?
- Identify a few or even just one good bright target(s)
- Adapt orbit to continuously view targets continuously or for ~month-long campaigns
- Just point and stare and collect tens of thousands of spectra
- Model the in- and out- of spectra transits

8. What are the best targets?
- Would need to be bright, as already shown
- Would need to have planets with close-in orbits to increase the number of transits viewable in a campaign duration
- Such planets would be highly irradiated.

9. What is the unique, new science
- Measuring the near-IR spectrum of a planet (still few examples)
- Detect water from molecular bands in the near-IR

10. Why is transit spectroscopy so hard?
- The integrated solid angle of the planet atmosphere is generally much less than the star's disk, yielding low signal-to-noise ratio of the atmospheric signal.
- Many consecutive observations are required to build up enough signal to noise ratio
- Systematic changes in the spectrograph can confound the disentangling of minuscule instrumental artifacts and planet-transit signals
- Unmodelled astrophysical variations (starspots / plages) also confound retrieval of planetary atmosphere signals.

11. What will you actually do for APRA 2018?  Why are technical barriers important?
- We will develop to two key enabling technologies for a compact near-IR spectrograph cubesat concept.
- Bonded cross-dispersed silicon immersed grisms for high efficiency, high resolution, high bandwidth, compact spectrograph design.
- Extreme precision detector metrology to calibrate detector systematics, enabling us to foward model the data acquisition process.
- Silicon immersed gratings deliver higher spectral grasp for a given beam size than conventional optical materials.  The >2x size savings reduces cost for a given desired spectral resolution, making a compact, low-cost cubesat design possible.  However, bonded, silicon grisms have not yet been demonstrated in space (their single)

12. Who will do the work?  Where?

- MGS will develop the Si grism technology at either UTexas, MIT Lincoln Lab, NASA Ames, or JPL.
- JL will develop the extreme precision detector metrology at MIT Lincoln Lab
- MGS will supervise Ames personnel will build simulations of correlated instrumental noise impact on exoplanet atmospheric retrievals.

13. What are the costs?
- MGS 20% time in first year, 30% time in second year.
- JL XX% time
- Additional postdoc and students at 50-100% time for N~2-3 years
- Custom equipment for detector metrology: ~$50K
- Experimental IR detector FPA and custom ASIC and electronics: ~$250K ???
- Si boule acquisition, orienting, cutting, lab costs: ~$50K
- Total: ~$1.0 - 1.5M ??


14. What other things to say?
- Demonstrate the low-cost cubesat platform.
- Demonstrate resilient strategies for calibrating instruments in Earth orbit.
- Innovations in software and modeling.
- Simulate a transit spectrum retrieval
- See Caroline Morley simulations

15. Why is your team uniquely suited to carrying out the proposed research program?

- MGS led development of Si grisms and immersion gratings for his PhD.  He developed a Fabry-Perot inspired metrology technique for detecting gaps as small as 10 nm in bonded Si optics (cite Gully-Santiago; APRA NNX)

- JL has pioneered custom ASICs and advanced detector applications for over a decade.

16. What are intrapixel sensitivity variations, and how does they affect precision measurements?

- Although much emphasis has focused on detector *read noise*, less metrology effort has focused on intrapixel sensitivity variations (cite ongoing work at RIT / Dmitri).  These minute perturbations are the main confounding factor in the unbiased analysis of Kepler/K2 data, and will limit TESS similarly.  Other factors like detector rolling band, and cross-talk also hamper our unbiased analysis of Kepler data.  By characterizing these instrumental signals, we will be able to forward model the generating process of our data acquisition.  We will have excellent generative models of both instrumental systematics, and excellent data-driven models of the stellar atmosphere from many repeated observations of transit-free spectra.  The final ingredient, and the desired science target---transiting exoplanet atmospheric spectra---will be inferred with a powerful probabilistic graphical model capable of detecting weak, correlated signals from messy residuals (cite Czekala et al. 2015, 2017, Gully-Santiago et al. 2017.)
