/*
    This SQL query is to select stars with high-fidelity measurements from the
    Gaia-ESO Survey which can be used for the calibration and/or verification of
    the second Gaia data release.

    Requested by: Andreas J. Korn (Uppsala)

    Author: Andrew R. Casey (Cambridge)
*/


/* How many "high fidelity" stars if we ignore their DR5 stellar parameters? */
SELECT  COUNT(DISTINCT cname) cname
  FROM  spectra
 WHERE  setup IN ('U580', 'U520', 'HR21', 'HR10') 
   AND  snr > 50
   AND  vel <> 'NaN'
   AND  e_vel <> 'NaN'
   AND  ges_type NOT LIKE 'AR%';

    /*
        SNR > 30: 30857
        SNR > 50: 20608
    */


/*  Also require the temperature to be in some sensible range, for other metrics
    to look sensible, and only use results from WG10 or WG11. */

SELECT  COUNT(DISTINCT s.cname) cname 
  FROM  spectra AS s,
        wg_recommended_results AS r 
 WHERE  r.cname = s.cname 
   AND  r.cname <> 'ssssssss-sssssss'
   AND  (
            (s.setup IN ('U580', 'U520') AND s.snr > 30)
        OR  (s.setup IN ('HR10', 'HR21') AND s.snr > 50)
        )
   AND  r.wg IN (10, 11)
   AND  s.vel <> 'NaN'
   AND  s.e_vel <> 'NaN'
   AND  s.ges_type NOT LIKE 'AR%'
   AND  r.teff > 4000
   AND  r.teff < 6500 
   AND  r.e_teff < 250 
   AND  r.nn_nodes_teff >= 2;


/*  Retrieve the required information. */

SELECT  DISTINCT ON (s.cname) s.ra, s.dec, s.cname
  FROM  spectra AS s,
        wg_recommended_results AS r 
 WHERE  r.cname = s.cname 
   AND  r.cname <> 'ssssssss-sssssss'
   AND  (
            (s.setup IN ('U580', 'U520') AND s.snr > 30)
        OR  (s.setup IN ('HR10', 'HR21') AND s.snr > 50)
        )
   AND  r.wg IN (10, 11)
   AND  s.vel <> 'NaN'
   AND  s.e_vel <> 'NaN'
   AND  s.vel > -500 AND s.vel < 500
   AND  s.ges_type NOT LIKE 'AR%'
   AND  r.teff > 4000
   AND  r.teff < 6500
   AND  r.e_teff < 250 
   AND  r.nn_nodes_teff >= 2;