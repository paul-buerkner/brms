## Test environments
* local Windows 8.1 install, R 3.2.0, R 3.3.0
* win-builder (devel and release)
* ubuntu 12.04, R 3.2.0
* local OS X install, R 3.2.0

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
  
* checking CRAN incoming feasibility ... NOTE   
  Suggests or Enhances not in mainstream repositories: rstan  
  Availability using Additional_repositories specification:  
  rstan   yes   http://rstan.org/repo/
    
  To my information, package rstan is not fully available on CRAN, yet, 
  because of some little vestiges of C++11 in its code. This should
  usually not cause any problems, although it is not offically part of C++03.
  
## Downstream dependencies
There are currently no downstream dependencies for this package.