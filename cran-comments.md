# Submission notes
* Update to version 0.3.2
* Replaced deprecated arma::is_finite(val) with std::isfinite(val)

## Reverse Dependency
Checked with tools::dependsOnPkgs() and tools::check_packages_in_dir()
* no reverse dependencies

## Test environments
Checked and passed using win-builder.r-project.org:
* windows (R-devel)
* windows (R-release)
* windows (R-old release)

Check and passed using R-CMD-check github actions:
* macos-latest (release)
* ubuntu-latest (devel)
* ubuntu-latest (oldrel-1)
* ubuntu-latest (release)
* windows-latest (release)

