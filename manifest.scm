;;; Copyright © 2024 Samuel Ortion <samuel@ortion.fr>
;;; guix manifest to install all dependencies for WGDdetector.smk
;;;
;;; use guix package --manifest=manifest.scm

(specifications->manifest
 '(
   "snakemake"
   "python"
   "r"
   "perl"
   "gawk"
   "blast+"
   ))
   ; "mmseqs2"
   ; "fastme" ;; not available on guix apparently
