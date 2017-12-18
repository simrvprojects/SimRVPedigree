##Resubmission
In this resubmission I have:

 * added a bioRxiv reference in the 'Description' field of the DESCRIPTION file in the form:
 authors (year) <bioRxiv:...>
 
 Note: After adding the reference R CMD check results issues a note because the reference does not end in a period.

## Test environments
* local Windows OS install, R 3.4.3
* ubuntu 14.05.5 LTS (on travis-ci), R 3.4.2
* local OS X install, R 3.3.3 
  - Tested on a campus MAC, which unfortunately I did not have permission to update to more recent version of R.

## R CMD check results
0 errors | 0 warnings | 1 note

 * 1 note: 
 Malformed Description field: should contain one or more complete sentences.


## Downstream Dependencies
Currently, there are no downstream dependencies for this package.
