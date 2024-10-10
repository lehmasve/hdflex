# Submission notes
* Update to version 0.3.0

## Reverse Dependency
Checked with tools::dependsOnPkgs() and tools::check_packages_in_dir()
* no reverse dependencies

## Test environments
Checked and passed using rhub v2.0.0 in the following environments:
* R-hubv2 linux R-* (any version)
* R-hubv2 macos R-* (any version)
* R-hubv2 macos-arm64 R-* (any version)
* R-hubv2 windows R-* (any version)
* R-hubv2 atlas R-devel
* R-hubv2 clang-asan R-devel
* R-hubv2 clang19 R-devel
* R-hubv2 ubuntu-clang R-devel
* R-hubv2 ubuntu-gcc12 R-devel
* R-hubv2 ubuntu-next R-4.4.1
* R-hubv2 ubuntu-release R-4.4.1

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

## Notes
* C++11 specification needed for parallelization over RcppThread







and 