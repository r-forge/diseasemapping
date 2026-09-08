## CRAN submission notes (draft)

Tarballs (ready to submit after your review / SVN commit / r-universe green):

* `pkg/geostatsp_2.1.1.tar.gz`
* `pkg/mapmisc_2.1.4.tar.gz`

### geostatsp 2.1.1 — local `R CMD check --as-cran`

Status: **1 NOTE** (expected)

* Suggests/Enhances `RandomFields` and `INLA` are not on CRAN;
  `INLA` is available via `Additional_repositories`.

Fixes relative to CRAN 2.0.11 / r-universe 2.1.0:

* Remaining `structure(.Dim=)` deprecation NOTE
* `inlaAvailable()` / `inlaSetThreads()` so broken INLA installs do not crash examples/tests
* `glgm()` returns a structured failure object (`parameters=NULL`) when INLA fails
* Test PDFs to `tempdir()`; packaging hygiene; terra documentation wording

### mapmisc 2.1.4 — local `R CMD check --as-cran`

Status: **OK**

* terra docs (including Nominatim `geocode`)
* HTTPS tile URLs + identifying User-Agent; dead providers pruned
* API-key note for Stadia/Thunderforest in `openmap` docs
* Packaging hygiene; `geocodeOld` removed; historical `inst/extsrc` marked

### Next steps (for you)

1. Review `svn status` under `pkg/geostatsp` and `pkg/mapmiscTerra`
2. Commit to R-Forge SVN when satisfied
3. Watch https://eborgnine.r-universe.dev/builds for green builds
   (especially linux-arm64 and macOS-oldrel for geostatsp)
4. Submit tarballs to CRAN
