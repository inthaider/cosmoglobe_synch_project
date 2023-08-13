md_jibran_v00.dat:
     Data obtained from /mn/stornext/d16/cmbco/AST9240/2022/data/md_hke_tmp.dat used for the collective run.

md_jibran_v01.dat:
     0.4-Haslam:
          std_mp        --- dec to 0.000E+00 to keep fixed
          std_dp        --- dec to 0.100E_08 to make closer to 0 ???
     0.080-LWA1:
          prior_mean_*  --- values obtained from healpy md fitting after doing approx galactic cut
          std_mp        --- inc to 0.100E+04 as suggested by Daniel (to see change in sampling)
          std_dp        --- inc to 0.100E+03 as suggested by Daniel (to see change in sampling)


md_jibran_v02.dat:
     0.080-LWA1:
          prior_mean_*  --- values obtained from healpy md fitting after doing approx galactic cut
          std_mp        --- dec back to 0.100E-02
          std_dp        --- dec back to 0.100E-05