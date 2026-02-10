J/MNRAS/402/620     SDSS WD main-sequence binaries  (Rebassa-Mansergas+, 2010)
================================================================================
Post-common envelope binaries from SDSS.
VII. A catalogue of white dwarf-main sequence binaries.
    Rebassa-Mansergas A., Gansicke B.T., Schreiber M.R., Koester D.,
    Rodriguez-Gil P.
   <Mon. Not. R. Astron. Soc., 402, 620-640 (2010)>
   =2010MNRAS.402..620R
================================================================================
ADC_Keywords: Stars, double and multiple ; Stars, white dwarf ; Stars, dwarfs ;
              Radial velocities ; Stars, masses
Keywords: stars: AGB and post-AGB - binaries: close - binaries: spectroscopic -
          stars: low-mass, brown dwarfs - white dwarfs

Abstract:
    We present a catalogue of 1602 white-dwarf-main-sequence (WDMS)
    binaries from the spectroscopic Sloan Digital Sky Survey Data Release
    6 (SDSS DR6). Among these, we identify 440 as new WDMS binaries. We
    select WDMS binary candidates by template fitting all 1.27 million DR6
    spectra, using combined constraints in both chi^2^ and signal-to-noise
    ratio. In addition, we use Galaxy Evolution Explorer (GALEX) and UKIRT
    Infrared Sky Survey (UKIDSS) magnitudes to search for objects in which
    one of the two components dominates the SDSS spectrum. We use a
    decomposition/fitting technique to measure the effective temperatures,
    surface gravities, masses and distances to the white dwarfs, as well
    as the spectral types and distances to the companions in our catalogue

Description:
    We have developed a procedure based on chi^2^ template fitting in
    order to automatically identify WDMS binary candidates from the SDSS
    DR6 spectroscopic data base.

File Summary:
--------------------------------------------------------------------------------
 FileName  Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe        80        .   This file
table2.dat    66      204   Updated classification of the 204 objects from
                             Silvestri et al. (2007, Cat <J/AJ/134/741>)
                             which are not considered as WDMS binaries by us
table3.dat    61       81   Updated classification of the 81 objects from
                             Heller et al. (2009, Cat. <J/A+A/496/191>)
                             which are not considered as WDMS binaries by us
table5.dat   242     1602   The complete catalogue. Coordinates and
                             GALEX-SDSS-UKIDSS magnitudes for the 1602 WDMS
                             binaries and candidates are also included
table7.dat   118     3592  *White-dwarf masses, effective temperatures,
                             surface gravities, spectral types and distances
                             of the 1602 WDMS binaries in our catalogue, as
                             determined from spectral modelling
table9.dat    68     1225   Radial velocities measured from the NaI
                             {lambda}{lambda}8183.27,8194.81 doublet and the
                             H{alpha} emission for 1068 systems in our catalogue
--------------------------------------------------------------------------------
Note on table7.dat: We list both the hot and the cold solutions, with the
     preferred solution given in the first line for each spectrum. The
     other solution is given for completeness.
--------------------------------------------------------------------------------

See also:
 J/MNRAS/382/1377 : SDSS WD main-sequence binaries (Rebassa-Mansergas+, 2007)
 J/ApJS/167/40    : SDSS4 confirmed white dwarfs catalog (Eisenstein+, 2006)
 J/AJ/134/741     : New close binary systems from SDSS DR5 (Silvestri+, 2007)
 J/A+A/486/843    : White dwarf-red dwarf binaries in SDSS (Augusteijn+, 2008)
 J/A+A/496/191    : Spectral analysis of 636 SDSS WD-M binaries (Heller+, 2009)

Byte-by-byte Description of file:table[23].dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  4  A4    ---     ---       [SDSS]
   5- 23  A19   ---     SDSS      SDSS designation (JHHMMSS.ss+DDMMSS.s)
  25- 66  A42   ---     Com       Comments
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table5.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  4  A4    ---     ---       [SDSS]
   5- 23  A19   ---     SDSS      SDSS designation (JHHMMSS.ss+DDMMSS.s)
  25- 34  F10.6 deg     RAdeg     Right ascension in decimal degrees (J2000)
  36- 45  F10.6 deg     DEdeg     Declination in decimal degrees (J2000)
  47- 54  F8.5  mag     NUV       ?=0 GALEX NUV magnitude (177-300nm)
  56- 62  F7.5  mag   e_NUV       ?=0 rms uncertainty on NUV
  64- 71  F8.5  mag     FUV       ?=0 GALEX FUV magnitude (135-175nm)
  73- 79  F7.5  mag   e_FUV       ?=0 rms uncertainty on FUV
  81- 88  F8.5  mag     umag      ?=0 SDSS u magnitude
  90- 96  F7.5  mag   e_umag      ?=0 rms uncertainty on umag
  98-105  F8.5  mag     gmag      ?=0 SDSS g magnitude
 107-113  F7.5  mag   e_gmag      ?=0 rms uncertainty on gmag
 115-122  F8.5  mag     rmag      ?=0 SDSS r magnitude
 124-130  F7.5  mag   e_rmag      ?=0 rms uncertainty on rmag
 132-139  F8.5  mag     imag      ?=0 SDSS i magnitude
 141-148  F8.5  mag   e_imag      ?=0 rms uncertainty on imag
 150-157  F8.5  mag     zmag      ?=0 SDSS z magnitude (1)
 159-166  F8.5  mag   e_zmag      ?=0 rms uncertainty on zmag
 168-175  F8.5  mag     Ymag      ?=0 UKIDSS Y magnitude (1)
 177-184  F8.6  mag   e_Ymag      ?=0 rms Yncertainty on Ymag
 186-193  F8.5  mag     Jmag      ?=0 UKIDSS J magnitude (1)
 195-201  F7.5  mag   e_Jmag      ?=0 rms uncertainty on Jmag
 203-210  F8.5  mag     Hmag      ?=0 UKIDSS H magnitude (1)
 212-218  F7.5  mag   e_Hmag      ?=0 rms uncertainty on Hmag
 220-227  F8.5  mag     Kmag      ?=0 UKIDSS K magnitude (1)
 229-235  F7.5  mag   e_Kmag      ?=0 rms uncertainty on Kmag
 237-242  A6    ---     Type      Object type
--------------------------------------------------------------------------------
Note (1): UKIDSS filters 50% cut-off wavelengths ({mu}m):
      -------------------
       Z    0.83   0.925
       Y    0.97   1.07
       J    1.17   1.33
       H    1.49   1.78
       K    2.03   2.37
      -------------------
     (details at http://www.ukidss.org/technical/instrument/filters.html)
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table7.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label    Explanations
--------------------------------------------------------------------------------
       1  I1    ---       Sol      [0/1]? Adopted (1) or ruled-out (0) solution
   3-  6  A4    ---       ---      [SDSS]
   7- 25  A19   ---       SDSS     SDSS designation (JHHMMSS.ss+DDMMSS.s)
  27- 32  A6    ---       Type     Object type
  34- 38  I5    ---       MJD      SDSS MJD identifier
  40- 43  I4    ---       Plate    SDSS Plate identifier
  45- 47  I3    ---       Fiber    SDSS Fiber identifier
  49- 53  I5    K         Teff     ?=0 Effective temperature
  55- 59  I5    K       e_Teff     ?=0 rms uncertainty on Teff
  61- 64  F4.2  [cm/s2]   logg     ?=0 Surface gravity
  66- 69  F4.2  [cm/s2] e_logg     ?=0 rms uncertainty on logg
  71- 75  F5.2  solMass   Mwd      ?=0.00 Mass of the white dwarf
  77- 80  F4.2  solMass e_Mwd      ?=0.00 rms uncertainty on Mwd
  82- 85  I4    pc        dwd      ?=0 Distance of the white dwarf
  87- 90  I4    pc      e_dwd      ?=0 rms uncertainty on dwd
  92- 93  I2    --        MType    ?=-1 Sub-type of the M dwarf
  95- 98  I4    pc        dsec     ?=0 Distance of the secondary
 100-103  I4    pc      e_dsec     ?=0 rms uncertainty on dsec
 105-113  A9    ---       Flag     Flag(s) for previously studied systems and
                                    resolved WDMS binary pairs (2)
 115-118  A4    ---       Notes    ? Individual notes (3)
--------------------------------------------------------------------------------
Note (2): Flags as follows:
      e = Eisenstein et al. (2006, Cat. <J/ApJS/167/40>)
      s = Silvestri et al. (2007, Cat <J/AJ/134/741>)
      a = Augusteijn et al. (2008, Cat <J/A+A/486/843>)
      h = Heller et al. (2009, Cat. <J/A+A/496/191>)
     re = resolved WDMS binary pairs
Note (3): Notes as follows:
      1 = Nearby Galaxy
      2 = GALEX magnitudes changes the solution obtained from the fit
      3 = No SDSS image available
      4 = Star/s nearby
      5 = PCEB
      6 = One of the spectra comes from a SEGUE plate
      7 = Irradiated secondary star
      8 = Broken spectrum
      9 = Star superimposed in the SDSS image,triple?
     10 = Wrong flux calibration in at least one spectrum
     11 = Cataclysmic Variable?
     12 = Wide binary?
     13 = Peculiar secondary star
     14 = Detached magnetic Cataclysmic Variable
     15 = Bad SDSS image
     16 = Low accretion rate polar
     17 = PG1159
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table9.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label   Explanations
--------------------------------------------------------------------------------
   1-  4  A4    ---     ---     [SDSS]
   5- 23  A19   ---     SDSS    SDSS designation (JHHMMSS.ss+DDMMSS.s)
      24  A1    ---     ---     [_]
      25  I1    ---   m_SDSS    [1/5]? Multiplicity index on SDSS
  27- 38  F12.4 d       HJD     Heliocentric Julian date
  40- 46  F7.2  km/s    RVNa    ?=0 Radial velocities measured from the
                                     NaI {lambda}{lambda}8183.27,8194.81 doublet
  48- 52  F5.2  km/s  e_RVNa    ?=0 rms uncertainty on RVNa
  54- 60  F7.2  km/s    RVHa    ?=0 Radial velocities measured from the H{alpha}
                                     emission
  62- 66  F5.2  km/s  e_RVHa    ?=0 rms uncertainty on RVHa
      68  A1    ---     Rel     [y/n] Radial velocity values obtained from
                                 spectra that are (y) and that are not (n)
                                 combined from individual exposures taken on
                                 different nights
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

References:
    Rebassa-Mansergas et al., Paper I    2007MNRAS.382.1377R, <J/MNRAS/382/1377>
    Schreiber et al,          Paper II   2008A&A...484..441S
    Rebassa-Mansergas et al., Paper III  2008MNRAS.390.1635R
    Nebot Gomez-Moran et al., Paper IV   2009A&A...495..561N
    Pyrzas et al.,            Paper V    2009MNRAS.394..978P
    Schwope et al.,           Paper VI   2009A&A...500..867S
================================================================================
(End)                                      Patricia Vannier [CDS]    19-Mar-2010
