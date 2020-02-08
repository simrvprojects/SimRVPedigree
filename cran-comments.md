##Resubmission February 2020
- Updated the DESCRIPTION file:
    - Updated the version number to 0.4.3 to reflect:
         - new functionality that allows users to simulate diseases with multiple disease subtypes.
         - removal of previously deprecated argument birth_range from the sim_life, sim_ped, and sim_RVped functions.
         - fixed all improper uses of class() == *, and instead switched to using inherits() or is().
         Many thanks to Kurt Hornik for again pointing out this deficiency. 

Test Environments prior to resubmission
* https://win-builder.r-project.org/
* local Windows OS install, R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0


##Resubmission October 2018
In this resubmission I have:
- Fixed a bug in the find_mrca function and improved the tests for this function. 

- We have improved the process used to simulate reproduction.  The new process better models age-specific fertility. As a result, the birth_range argument has been deprecated in the sim_life, sim_ped, and sim_RVped functions.  To prepare users for this change, in this resubmission, users are issued a warning if they alter the default NULL value of this argument.

- Updated the DESCRIPTION file:
    - Updated the version number to 0.3.0 to reflect:
         - the updated reproduction process, which better models age-specific fertility.
         - the deprecated argument birth_range from the sim_life, sim_ped, and sim_RVped functions.
         - the bug fix in the find_mrca function,
    - Reverted the Depends field to R (>= 3.4.0) to fix error on package check page for flavor r-oldrel-osx-x86_64.

Test Environments prior to resubmission
* local Windows OS install, R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0


##Resubmission August 2018
In this resubmission I have:
- Fixed a bug in the reassign_gen function and improved the tests for this function. 

- Added a new argument, called first_diagnosis, to the sim_RVped function.  This argument allows users to implement new ascertainment criteria. 
  NOTE: the default setting of the new argument ensures backwards compatibility. 

- Added a function called find_mrca, which finds the most recent common ancestor of two pedigree members.


- Updated the DESCRIPTION file:
    - Updated the Depends field to R (>= 3.5.0),
    - Updated the version number to 0.2.0 to reflect:
         - the bug fix in the reassign_gen function,
         - the addition of the new argument in the sim_RVped function, called first_diagnosis, and
         - the addition of a new function called reduce_to_affected.

Test Environments prior to resubmission
* local Windows OS install, R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0


##Resubmission
In this resubmission I have:

 * fixed the reference in the 'Description' field of the DESCRIPTION file to be of the form:
 authors (year) <doi...>, and

 * changed the 'Depends' field in the DESCRIPTION file to R (>= 3.4.0). Previously, it was R (>= 3.4.2).

##Resubmission
In this resubmission I have updated a reference in AgeSpecififc_Hazards.Rd, which had prompted the following error message:
 
 Found the following (possibly) invalid URLs:
   URL: https://www.ssa.gov/oact/NOTES/as120/TOC.html
     From: man/AgeSpecific_Hazards.Rd
           inst/doc/SimRVPedigree.html
     Status: Error
     Message: libcurl error code 35:
               Unknown SSL protocol error in connection to www.ssa.gov:443
               
 * The updated reference does not contain a URL, and should therefore resolve this error.

##Resubmission
In this resubmission I have:

 * fixed the reference in the 'Description' field of the DESCRIPTION file to be of the form:
 authors (year) <https...>
 
## Test environments
* local Windows OS install, R 3.4.3
* ubuntu 14.05.5 LTS (on travis-ci), R 3.4.2
* local OS X install, R 3.3.3 
  - Tested on a campus MAC, which unfortunately I did not have permission to update to more recent version of R.

## R CMD check results
0 errors | 0 warnings | 1 note

Possibly mis-spelled words in DESCRIPTION:
  Jinko (12:9)
  Nieuwoudt (11:15)
  
These are not misspelled words, these are names of authors. 

## Downstream Dependencies
Currently, there are no downstream dependencies for this package.
