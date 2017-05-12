KPL/FK

Frame (FK) SPICE kernel file for JUICE-specific generic frames
===============================================================================

   This frames kernel defines a number of generic frames used by JUICE
   mission for Science opportunities identification, data analysis and
   scientific research. These frames are currently not ``built'' into
   the SPICE toolkit.
   
   These frames are sorted in two groups: those that are JUICE mission 
   specific and those that are Jupiter system generic. The first group
   contains the frames defined by and for the JUICE mission, while the
   second provides the frames that are commonly accepted by the scientific
   community for the Jupiter system and its satellites.
   
   The IAU body-fixed rotational frames for Jupiter, Callisto, Europa and
   Ganymede is an exception to this grouping, as they are provided in a
   separate PCK kernel file. 
   

Version and Date
-------------------------------------------------------------------------------

   Version 0.3 -- July 08, 2016 -- Marc Costa Sitja, ESAC/ESA

      The updates of this version have been motivated by the feedback
      of the JUICE Working Group 3 of the Science Operations Working
      Group iterated with NAIF.

      Updated reference [14] and added references [15], [16] and [17].
      
      Added JUICE_JUPITER_IF_J2000 frame as described in reference [15] 
      and as provided by Boris Semenov.
      
      Added a transformed version of CALLISTO_JUPITER_ORB, EUROPA_JUPITER_ORB  
      and GANYMEDE_JUPITER_ORB frames to use the same convention as the
      community that models the local moon-magnetosphere interaction. The
      transformation applied is x'=-y, y'=x, z'=z. The new frames are  
      JUICE_CALLISTO_PHI_ORB, JUICE_EUROPA_PHI_ORB and JUICE_GANYMEDE_PHI_ORB.

      Updated JUICE_JUPITER_MAG frame with dipole and prime meridians
      definitions from reference [16] as used by the JUNO mission.
      Changed name to JUICE_JUPITER_MAG_SIII.

      Added JUICE_JUPITER_SIII_RH frame as described in references 
      [14] and [17].

   Version 0.2 -- June 04, 2016 -- Marc Costa Sitja, ESAC/ESA

      Updated all NAIF ID codes from -907* to -28* since the JUICE 
      spacecraft NAIF ID has been updated from -907 to -28. 

      Added comments for future updates.

   Version 0.1 -- April 26, 2016 -- Jorge Diaz del Rio (ODC Space)
   
      Modified frames' names to make them mission specific by adding JUICE_*
      Changed secondary axis in Jupiter Magnetic frame to be X instead of Y.
      Minor edit to comments and section descriptions.

   Version 0.0 -- February 5, 2016 -- Jorge Diaz del Rio (ODC Space)
   
      Initial version.

      
References
-------------------------------------------------------------------------------

   [1]   "Frames Required Reading"
   
   [2]   "Kernel Pool Required Reading"
   
   [3]   Kivelson, M., Russell, C., 1995. "Introduction to Space
         Physics." Cambridge Univ. Press, Cambridge,UK.

   [4]   http://stereo.sr.unh.edu/data/PLASTIC_Resources/
         stereo_coordinates.pdf

   [5]   Weiduo Hu, "Fundamental Spacecraft Dynamics and Control,"
         Wiley 2015

   [6]   Bentley, R.D., Hapgood, M. and Thompson, W.T., "Coordinate
         Systems," 22-Mar-2010, HELIO-UCL=N3-003-TN.
         http://www.helio-vo.eu/documents/public/
         HELIO_Coordinates_100322.pdf

   [7]   Franz and Harper. (2002) "Heliospheric Coordinate Systems,"
         Space Science, 50, 217ff.

   [8]   Burlaga, L. F., 1984. "MHD processes in the outer
         heliosphere." Space Sci. Rev. 39, 255-316.

   [9]   "Space Station Reference Coordinate Systems," Revision F,
         26-Oct-2001

   [10]  Li, H., "Geostationary Satellites Collocation," Springer, 2014

   [11]  Russell, C.T., "Geophysical Coordinate Transformations,"
         Cosmic Electrodynamics, 2, 184-196, 1971
         
   [12]  Bhavnani, K.H. and Vancour, R.P. (1991), "Coordinate Systems for
         Space and Geophysical Applications"
         
   [13]  NASA/Marshall Solar Physics: 
         http://solarscience.msfc.nasa.gov/SolarWind.shtml
         
   [14]  Seidelmann, P.K., Divine, N. "Evaluation of Jupiter longitudes in 
         System III (1965)", Geophyisical Research Letters, Volume 4, 
         Issue 2, Pages 65–68, 1977

   [15]  Liu X., Sachse, M., Spahn, F., Schmidt J., "Dynamics and 
         distribution of Jovian dust ejected from the Galilean satellites",
         doi: 10.1002/2016JE004999, American Geophysical Union, 2016

   [16]  Connerneyand, J.E.P., Acuna, M.H., Ness N.F., Satoh, T.,  
         "New models of Jupiter's magnetic field constrained by the Io 
         flux tube footprint", Journal of Geophysical Research, Vol. 103, 
         No. A6, 929-939, 1998

   [17]  Bagenal F., Wilson, R.J., "Jupiter Coordinate Systems", 
         http://lasp.colorado.edu/home/mop/missions/juno/coordinates/,
         accessed 24/06/16.

      
Contact Information
-------------------------------------------------------------------------------

   If you have any questions regarding this file contact SPICE support at
   ESAC:

           Marc Costa Sitja
           (+34) 91-8131-457
           mcosta@sciops.esa.int, esa_spice@sciops.esa.int
           
   or NAIF at JPL:
   
           Boris Semenov
           (818) 354-8136
           Boris.Semenov@jpl.nasa.gov
           

Implementation Notes
-------------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make
   use of this frame kernel must "load" the kernel normally during 
   program initialization. Loading the kernel associates the data items 
   with their names in a data structure called the "kernel pool". The 
   routine that loads a kernel into the pool is shown below:
                                                                               
      FORTRAN: (SPICELIB)

         CALL FURNSH ( frame_kernel_name )

      C: (CSPICE)

         furnsh_c ( frame_kernel_name );

      IDL: (ICY)

         cspice_furnsh, frame_kernel_name
         
      MATLAB: (MICE)
      
         cspice_furnsh ( 'frame_kernel_name' )

   This file was created and may be updated with a text editor or word
   processor

   
SPICE Frame names and NAIF ID Codes
-------------------------------------------------------------------------------
 
   The following generic frames are defined in this kernel file:

      SPICE Frame Name              Long-name
      --------------------------    -------------------------------------------
        
   JUICE mission specific generic frames(*):
  
      JUICE_JUPITER_IF_J2000        Jupiter Inertial Frame at J2000
      JUICE_JUPITER_SIII_RH         Jupiter System III Right-Handed    
  
      JUICE_JUPITER_MAG_SIII        Jupiter Magnetic System III Right-Handed
      JUICE_JUPITER_BEQXS           Jovian-centric Equatorial Solar
      JUICE_JUPITER_BSM             Jupiter Solar Magnetospheric
      JUICE_JUPITER_SJC             Solar Joviancentric Coordinates
      JUICE_JSM                     Jupiter Solar Magnetic
      JUICE_JSW                     Jupiter Solar Wind
      JUICE_JSWM                    Jupiter Solar Wind Magnetospheric
  
      JUICE_CALLISTO_PHI_ORB        Callisto Orbital Transformed
      JUICE_EUROPA_PHI_ORB          Europa Orbital Transformed
      JUICE_GANYMEDE_PHI_ORB        Ganymede Orbital Transformed
  
      JUICE_HGRTN                   JUICE Heliocentric Radial-Tangential-Normal
      JUICE_SUN_RTN                 Sun-orbit JUICE Radial-Tangential-Normal
      JUICE_JUPITER_DM              Jupiter JUICE Dipole Meridian
  
   Jupiter system generic frames(**):
  
      JUPITER_MEQUD                 Jupiter Mean Equator of date
      JUPITER_SUN_EQU               Jupiter Solar Equatorial
      JUPITER_SUN_ORB               Jupiter Solar Orbital
        
      JUPITER_CALLISTO_BCSF         Jupiter-centric Callisto-following    
      JUPITER_EUROPA_BCSF           Jupiter-centric Europa-following      
      JUPITER_GANYMEDE_BCSF         Jupiter-centric Ganymede-following    
  
      CALLISTO_JUPITER_ORB          Callisto Orbital
      EUROPA_JUPITER_ORB            Europa Orbital
      GANYMEDE_JUPITER_ORB          Ganymede Orbital 


   (*)  The definition of the dipole direction of these frames needs to
        be discussed within the Science Operations Working Team Working
        Group 3 in coordination with Juno and Europa Multi Fly-by Mission.
        The name for JUICE_SUN_RTN also needs to be discussed.

   (**) These frames are commonly used by other missions for data analysis
        and scientific research. In the future NAIF may include them
        in their official generic frames kernel for the Jupiter system.
        When this happens the frames will be removed from this kernel.
	      

   These frames have the following centers, frame class and NAIF
   IDs:
   
      SPICE Frame Name              Center                 Class    NAIF ID
      ---------------------------   ---------------------  -------  ----------
      JUICE_JUPITER_IF_J2000        JUPITER                FIXED        -28970
      JUICE_JUPITER_SIII_RH         JUPITER                FIXED        -28971
  
      JUICE_JUPITER_MAG_SIII        JUPITER                FIXED        -28972
      JUICE_JUPITER_BEQXS           JUPITER                DYNAMIC      -28973
      JUICE_JUPITER_BSM             JUPITER                DYNAMIC      -28974
      JUICE_JUPITER_SJC             JUPITER                DYNAMIC      -28975
      JUICE_JSM                     JUPITER                DYNAMIC      -28976
      JUICE_JSW                     JUPITER                FIXED        -28977
      JUICE_JSWM                    JUPITER                FIXED        -28978
  
      JUICE_CALLISTO_PHI_ORB        JUPITER                DYNAMIC      -28980
      JUICE_EUROPA_PHI_ORB          JUPITER                DYNAMIC      -28981
      JUICE_GANYMEDE_PHI_ORB        JUPITER                DYNAMIC      -28982
  
      JUICE_HGRTN                   SUN                    DYNAMIC      -28990
      JUICE_SUN_RTN                 JUICE                  DYNAMIC      -28991
      JUICE_JUPITER_DM              JUPITER                DYNAMIC      -28992
  
      JUPITER_MEQUD                 JUPITER                DYNAMIC   500599000
      JUPITER_SUN_EQU               JUPITER                DYNAMIC   500599001
      JUPITER_SUN_ORB               JUPITER                DYNAMIC   500599002
        
      JUPITER_CALLISTO_BCSF         JUPITER                DYNAMIC   500599006
      JUPITER_EUROPA_BCSF           JUPITER                DYNAMIC   500599004
      JUPITER_GANYMEDE_BCSF         JUPITER                DYNAMIC   500599005
     
      CALLISTO_JUPITER_ORB          CALLISTO               DYNAMIC   500504000
      EUROPA_JUPITER_ORB            EUROPA                 DYNAMIC   500502000
      GANYMEDE_JUPITER_ORB          GANYMEDE               DYNAMIC   500503000


   The keywords implementing these frame definitions are located in the 
   "JUICE Mission Specific Generic Frame Definitions" and "JUPITER System
   Generic Frame Definitions" sections.
               
      
General Notes About This File
-------------------------------------------------------------------------------

   About Required Data:
   --------------------
   Most of the dynamic frames defined in this file require at least one
   of the following kernels to be loaded prior to their evaluation, 
   normally during program initialization:

     - Planetary ephemeris data (SPK), i.e. de403, de405, etc;
     - Planetary constants data (PCK);
     - Earth generic frames definitions (FK).

   Note that loading different kernels will lead to different
   orientations of the same frame at a given epoch, providing different
   results from each other, in terms of state vectors referred to these 
   frames.
 
 
   About Implementation:
   ---------------------
   The SPICE frames defined within this file and their corresponding 
   references in literature might not be equivalent, both due to 
   variations in the SPICE kernels on which the SPICE frame depends, 
   and due to possible differences in both the frame's definition and 
   implementation (e.g. GSE can be defined using the instantaneous 
   orbital plane or mean ecliptic; the mean ecliptic is a function of 
   the ecliptic model). Please refer to each applicable frame 
   description section for particular details on the current SPICE
   kernel implementation.
   
 
JUICE Mission Specific Scientific Frame Definitions
-------------------------------------------------------------------------------

   This section contains the definition of the JUICE mission specific
   scientific frames.


Jupiter Intertial Frame at J2000 (JUICE_JUPITER_IF_J2000)
------------------------------------------------------------------------

   Definition:
   -----------   
   The Jupiter Intertial Frame at J2000 is defined as follows:

      -  +Z axis is parallel to Jupiter rotation axis 
         at J2000, pointing toward the North side of the 
         invariable plane;

      -  +X axis is aligned with the ascending node of the Jovian 
         orbital plane with the Jovian equator plane at J2000;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of 
         Jupiter.

   All vectors are geometric: no corrections are used.


   Remarks:
   --------
   This frame is defined as a fixed offset frame using constant vectors 
   as the specification method. The fixed offset for these vectors were 
   based on the following directions (that also define a two-vector 
   frame):

     - +Z axis along Right Ascension of 268.05720404270755 degrees and 
       Declination of 64.495809953395721 of Jupiter pole at J2000 in 
       J2000;

     - +X axis along Right Ascension of 139.90957808624046 degrees and 
       Declination of 16.419104754306645 with respect to J2000 at J2000,
       which corresponds to the RA/DEC of Jupiter instanteneous orbital 
       plane ascending node on Jupiter equator at J2000 in J2000; 

   This frame has been defined based on the IAU_JUPITER frame, whose 
   evaluation was based on the data included in the loaded PCK file.
   When evaluated the IAU 2009 constaints were used.

   In addition jup310.bsp has been used to compute the Jupiter 
   instanteneous orbital plane ascending node on Jupiter equator at 
   J2000 in J2000.


  \begindata

     FRAME_JUPITER_IF_J2000         = -28970
     FRAME_-28970_NAME              = 'JUPITER_IF_J2000'
     FRAME_-28970_CLASS             = 4
     FRAME_-28970_CLASS_ID          = -28970
     FRAME_-28970_CENTER            = 599

     TKFRAME_-28970_SPEC            = 'ANGLES'
     TKFRAME_-28970_RELATIVE        = 'J2000'
     TKFRAME_-28970_ANGLES          = ( 
                                         -88.057204042707568   
                                          25.504190046604286   
                                         -48.968728165793905
                                      )
     TKFRAME_-28970_AXES            = ( 
                                          3,
                                          2,   
                                          3
                                      )
     TKFRAME_-28970_UNITS           = 'DEGREES'

  \begintext


Jupiter System III Right-handed (JUICE_JUPITER_SIII_RH)
------------------------------------------------------------------------
   
   SPICE frame name and literature references
   ------------------------------------------
   Within the SPICE system, the Jupiter System III Right-handed is 
   referred as JUPITER_SIII_RH. In literature in the documentation 
   provided by the JUNO team (see [17]), this frame is referred as S3RH. 
   Please note that in the JUNO S3RH colatitude is used instead of 
   latitude.

   
   Definition:
   -----------
   The Jupiter System III frame is defined as follows (from [14]):

      -  +Z axis is parallel to Jupiter rotation axis, pointing 
         toward the North side of the invariable plane;
 
      -  +X axis is aligned with the ascending node of the Jovian 
         orbital plane with the Jovian equator plane;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Jupiter.


   Remarks:
   --------
   This frame is an alias for the IAU_JUPITER frame using a TK frame. 

   This frame is valid as long as the IAU_JUPITER is defined with 
   the PCK kernel pck00010.tpc. If the PCK kernel is updated and the
   Jupiter Pole and Prime Meridian definitions are updated this frame
   shall be defined by other means.
 

  \begindata
 
      FRAME_JUPITER_SIII_RH    =  -28971
      FRAME_-28971_NAME        = 'JUPITER_SIII_RH'
      FRAME_-28971_CLASS       =  4
      FRAME_-28971_CLASS_ID    =  -28971
      FRAME_-28971_CENTER      =  599
      
      OBJECT_599_FRAME         = 'JUPITER_SIII_RH'
 
      TKFRAME_-28971_RELATIVE  = 'IAU_JUPITER'
      TKFRAME_-28971_SPEC      = 'MATRIX'
      TKFRAME_-28971_MATRIX    = ( 1   0   0
                                   0   1   0
                                   0   0   1 )

  \begintext


Jupiter Magnetic System III Right-Handed (JUICE_JUPITER_MAG_SIII)
------------------------------------------------------------------------

   SPICE Frame Name and synonyms:
   ------------------------------
   Within the SPICE system, the Jupiter Magnetic System III Right-Handed 
   frame is referred as JUICE_JUPITER_MAG_SIII.


   Definition:
   -----------
   The Jupiter Magnetic System III Right-Handed frame is defined as 
   follows (from [17]):

    -  +Z axis is parallel to the magnetic dipole axis of epoch and aligned with
       the North magnetic pole of epoch;

    -  +X axis is defined with respect to the meridian where the magnetic and 
       geographic equators cross at longitude 69.2 degrees.

    -  +Y completes the right-handed frame;

    -  the origin of this frame is Jupiter's center of mass.
    
    
   Remarks:
   --------
   This frames is the Jupiter System III Right-handed frame but is tilted
   by the dipole approximation of the magnetic field of Jupiter (see [17]). 

   The north dipole location is the best-fit, centered dipole magnetic 
   moment vector M direction in the IAU_JUPITER frame (see [16]). M is
   defined to be titled 9.60 degrees away from the Jovian spin axis
   towards longitude 200.8 degrees West, therefore, the current values
   for the Jupiter north magnetic centered dipole planetocentric
   coordinates are:
      
      Longitude    =  159.2 (degrees)
      Latitude     =   80.5 (degrees)
 
   The magnetic longitude is defined with respect to the meridian where the 
   magnetic and geographic equators cross (where ThetaIII =0 deg and 
   ThehaRH = ThehaMAG =90 deg) at lIII=290.8 deg or lRH=69.2 deg.

  \begindata

      FRAME_JUICE_JUPITER_MAG_SIII        = -28972
      FRAME_-28972_NAME                   = 'JUICE_JUPITER_MAG_SIII'         
      FRAME_-28972_CLASS                  = 4
      FRAME_-28972_CLASS_ID               = -28972
      FRAME_-28972_CENTER                 = 599 
      TKFRAME_-28972_SPEC                 = 'ANGLES'
      TKFRAME_-28972_RELATIVE             = 'IAU_JUPITER'
      TKFRAME_-28972_ANGLES               = ( 0, -69.2, 9.5 )
      TKFRAME_-28972_AXES                 = ( 2,   3,   1   )
      TKFRAME_-28972_UNITS                = 'DEGREES'

  \begintext


Joviancentric Equatorial Solar frame (JUICE_JUPITER_BEQXS)
------------------------------------------------------------------------

   SPICE frame name, common names and other designators:
   -----------------------------------------------------
   Within the SPICE system, and for the JUICE mission, the Joviancentric
   Equatorial Solar frame is referred as JUICE_JUPITER_BEQXS. In the
   documentation provided by the JUICE MAG team (see [14]), it is referred as
   Joviancentric Solar Equatorial frame or JSEq.
   
   CAUTION: xSEQ coordinate systems are typically referenced to the planet Sun
   line and the solar equator. JUICE mission uses a lowercase q in the frame
   designator to distinguish this type of coordinate system, where the
   planetary equator is the plane of reference.
   

   Definition:
   -----------
   The Joviancentric Equatorial Solar frame is defined as follows (from [4]):

      -  X-Y plane is defined by the Jupiter equator of date. Therefore, the
         +Z axis is the primary vector and it is aligned to Jupiter's north
         pole of date;

      -  +X axis is the component of the Jupiter-Sun vector that is orthogonal
         to the +Z axis;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is Jupiter's center of mass.

   All the vectors are geometric: no aberration corrections are used.
   

   Uses and applications:
   ----------------------
   The Joviancentric Equatorial Solar frame is commonly used in the analysis
   of orbital trajectories where planetary rotation complicates the analysis
   but where knowledge of the equatorial plane is still required. Data can be
   analyzed for local time and north-south asymmetries in this reference frame.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different types of
   specifications for the primary and secondary vectors.

   The primary vector is defined as a constant vector in the IAU_JUPITER
   frame, which is a PCK-based frame. Therefore a PCK file containing the
   orientation constants for Jupiter has to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target position' vector.
   Therefore, the ephemeris data required to compute the Jupiter-Sun vector
   in J2000 frame have to be loaded before using this frame.
 
  
   Remarks:
   --------
   This frame is defined based on SPK data: different planetary ephemerides
   for Jupiter, Jupiter's Barycenter, the Sun and the Solar System Barycenter
   will lead to different frame orientation at a given time.
   
   This frame is also defined based on the IAU_JUPITER frame, whose evaluation
   is based on the data included in the loaded PCK file: different orientation
   constants for Jupiter's spin axis will lead to different frame orientation
   at a given time. 
   
   It is strongly recommended to indicate what data have been used in the
   evaluation of this frame when referring to it, i.e. JUICE_JUPITER_BEQXS using
   IAU 2009 constants and de403 ephemerides.
   
  \begindata

      FRAME_JUICE_JUPITER_BEQXS     =  -28973
      FRAME_-28973_NAME             = 'JUICE_JUPITER_BEQXS'
      FRAME_-28973_CLASS            =  5
      FRAME_-28973_CLASS_ID         =  -28973
      FRAME_-28973_CENTER           =  599
      FRAME_-28973_RELATIVE         = 'J2000'
      FRAME_-28973_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-28973_FAMILY           = 'TWO-VECTOR'
      FRAME_-28973_PRI_AXIS         = 'Z'
      FRAME_-28973_PRI_VECTOR_DEF   = 'CONSTANT'
      FRAME_-28973_PRI_FRAME        = 'IAU_JUPITER'
      FRAME_-28973_PRI_SPEC         = 'RECTANGULAR'
      FRAME_-28973_PRI_VECTOR       = ( 0, 0, 1 )
      FRAME_-28973_SEC_AXIS         = 'X'
      FRAME_-28973_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-28973_SEC_OBSERVER     = 'JUPITER'
      FRAME_-28973_SEC_TARGET       = 'SUN'
      FRAME_-28973_SEC_ABCORR       = 'NONE'

  \begintext


Jupiter Solar Magnetospheric (JUICE_JUPITER_BSM)
------------------------------------------------------------------------

   SPICE frame name and literature references
   ------------------------------------------
   Within the SPICE system, the Jupiter Solar Magnetospheric frame is 
   referred as JUICE_SUN_BSM. In the documentation provided by the
   JUICE MAG team (see [14]),, this frame is referred as JSM.


   Definition:
   -----------
   The Jupiter Solar Magnetospheric frame is defined, based on the
   definition in [11] for the Earth, as follows (from [14]):
   
      -  +X axis is the position of the Sun relative to Jupiter; 
         it's the primary vector and points from Jupiter to the 
         Sun;
         
      -  +Z axis is the projection of the magnetic centered
         dipole axis (positive North) of Jupiter onto the plane 
         perpendicular to the +X axis.
    
      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         the Jupiter.

   All vectors are geometric: no corrections are used.


   Uses and applications:
   ----------------------
   The Jupiter Solar Magnetospheric is used to analyze data affected
   by both the solar wind flow velocity and the rotating magnetic
   field of Jupiter, but where the solar wind forces still dominate
   (in the magnetosheath and near the magnetopause) (from [14]).
   
   
   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.
   
   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Sun vector in J2000 reference frame have to be loaded 
   before using this frame.
   
   The secondary vector is defined as a constant vector in the 
   JUICE_JUPITER_MAG body-fixed frame, which provides the orientation
   for the Jupiter North magnetic dipole and the Jupiter magnetic
   frame with respect to the IAU_JUPITER frame.   
   
   Remarks:
   --------
   This frame is defined based on SPK data: different planetary 
   ephemerides for Jupiter, Jupiter Barycenter, the Sun and the 
   Solar System Barycenter will lead to a different frame orientation 
   at a given time.
   
   The secondary vector of this frame is defined based on the 
   JUICE_JUPITER_MAG frame. Different implementations of this frame, i.e. 
   changes to the existing frame definition will lead to different
   orientations of the North geomagnetic centered dipole, and therefore
   to different projections of this vector on the frame's YZ plane,
   resulting on different orientations of this frame.
      
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, i.e. 
   JUICE_JUPITER_BSM using IAU 2009 constants and de405 ephemerides.
   
  \begindata

       FRAME_JUICE_JUPITER_BSM       =  -28974
       FRAME_-28974_NAME             = 'JUICE_JUPITER_BSM'
       FRAME_-28974_CLASS            = 5
       FRAME_-28974_CLASS_ID         = -28974
       FRAME_-28974_CENTER           =  599
       FRAME_-28974_RELATIVE         = 'J2000'
       FRAME_-28974_DEF_STYLE        = 'PARAMETERIZED'
       FRAME_-28974_FAMILY           = 'TWO-VECTOR'
       FRAME_-28974_PRI_AXIS         = 'X'
       FRAME_-28974_PRI_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
       FRAME_-28974_PRI_OBSERVER     = 'JUPITER'
       FRAME_-28974_PRI_TARGET       = 'SUN'
       FRAME_-28974_PRI_ABCORR       = 'NONE'
       FRAME_-28974_SEC_AXIS         = 'Z'
       FRAME_-28974_SEC_VECTOR_DEF   = 'CONSTANT'
       FRAME_-28974_SEC_SPEC         = 'RECTANGULAR'
       FRAME_-28974_SEC_VECTOR       = ( 0, 0, 1 )
       FRAME_-28974_SEC_FRAME        = 'JUICE_JUPITER_MAG_SIII'

  \begintext


Solar Joviancentric Coordinates frame (JUICE_JUPITER_SJC)
------------------------------------------------------------------------

   SPICE frame name, common names and other designators:
   -----------------------------------------------------
   Within the SPICE system, and for the JUICE mission, the Solar Joviancentric
   Coordinates frame is referred as JUICE_JUPITER_SJC. In the documentation
   provided by the JUICE MAG team (see [14]), it is referred as Solar
   Joviancentric Coordinates or SJ.
   
   
   Definition:
   -----------
   The Solar Joviancentric Coordinates frame is defined as follows (from [14]):
   
      -  +X axis is the position of the Sun relative to Jupiter; it's
         the primary vector and points from Jupiter to the Sun;

      -  +Z axis is the component of Jupiter's north pole of date 
         orthogonal to the +X axis;

      -  +Y axis completes the right-handed reference frame;

      -  the origin of this frame is Jupiter's center of mass.

   All the vectors are geometric: no aberration corrections are used.
   

   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different types
   of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' vector.
   Therefore, the ephemeris data required to compute the Jupiter-Sun vector
   in J2000 frame have to be loaded before using this frame.

   The secondary vector is defined as a constant vector in the IAU_JUPITER
   frame, which is a PCK-based frame. Therefore a PCK file containing
   the orientation constants for Jupiter has to be loaded before using
   this frame.


   Remarks:
   --------
   This frame is defined based on the IAU_JUPITER frame, whose evaluation is
   based on the data included in the loaded PCK file: different orientation
   constants for Jupiter's spin axis will lead to different frames.

   Since the primary vector of this frame is defined as an 'observer-target
   position' vector, the usage of different planetary ephemerides
   may lead to a different frame orientation at a given time.
   
   It is strongly recommended to indicate what data have been used in the
   evaluation of this frame when referring to it, i.e. JUICE_JUPITER_SJC using 
   IAU 2009 constants and de405 ephemerides.


  \begindata

      FRAME_JUICE_JUPITER_SJC       =  -28975
      FRAME_-28975_NAME             = 'JUICE_JUPITER_SJC'
      FRAME_-28975_CLASS            =  5
      FRAME_-28975_CLASS_ID         =  -28975
      FRAME_-28975_CENTER           =  599
      FRAME_-28975_RELATIVE         = 'J2000'
      FRAME_-28975_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-28975_FAMILY           = 'TWO-VECTOR'
      FRAME_-28975_PRI_AXIS         = 'X'
      FRAME_-28975_PRI_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-28975_PRI_OBSERVER     = 'JUPITER'
      FRAME_-28975_PRI_TARGET       = 'SUN'
      FRAME_-28975_PRI_ABCORR       = 'NONE'
      FRAME_-28975_SEC_AXIS         = 'Z'
      FRAME_-28975_SEC_VECTOR_DEF   = 'CONSTANT'
      FRAME_-28975_SEC_FRAME        = 'IAU_JUPITER'
      FRAME_-28975_SEC_SPEC         = 'RECTANGULAR'
      FRAME_-28975_SEC_VECTOR       = ( 0, 0, 1 )

  \begintext


Jupiter Solar Magnetic (JUICE_JSM)
------------------------------------------------------------------------

   Definition:
   -----------
   The Jupiter Solar Magnetic frame is defined, based on the definition
   provided in [11] for the Earth's solar magnetic frame, as follows:
   
      -  +Z axis is aligned with the North magnetic centered dipole of
         Jupiter;

      -  +X axis is the component of the Jupiter-Sun vector that is
         orthogonal to the +Z axis;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of Jupiter.

   All vectors are geometric: no corrections are used.
   
   
   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different types
   of specifications for the primary and secondary vectors.
   
   The primary vector is defined as a constant vector in the 
   JUICE_JUPITER_MAG body-fixed frame, which provides the orientation
   for the Jupiter North magnetic dipole and the Jupiter magnetic
   frame with respect to the IAU_JUPITER frame.   

   The secondary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Sun vector in J2000 reference frame have to be loaded
   before using this frame.


   Remarks:
   --------
   The primary vector of this frame is defined based on the 
   JUICE_JUPITER_MAG frame. Different implementations of this frame, i.e. 
   changes to the existing frame definition will lead to different
   orientations of the North geomagnetic centered dipole, but only when
   these data lead to different right ascension and declination
   coordinates of the North magnetic centered dipole vector of Jupiter.
        
   Since the secondary vector of this frame is defined as an 
   'observer-target position' vector, the usage of different planetary 
   ephemerides for Jupiter, Jupiter barycenter, the Sun and the Solar
   System barycenter may lead to different orientations of this frame, 
   but only when these data lead to different projections of the 
   Jupiter-Sun vector on the magnetic equator of Jupiter.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, i.e. 
   JUICE_JSM using IAU 2009 constants and de405 ephemerides.
   
  \begindata
  
       FRAME_JUICE_JSM               =  -28976
       FRAME_-28976_NAME             = 'JUICE_JSM'
       FRAME_-28976_CLASS            =  5
       FRAME_-28976_CLASS_ID         =  -28976
       FRAME_-28976_CENTER           =  599
       FRAME_-28976_RELATIVE         = 'J2000'
       FRAME_-28976_DEF_STYLE        = 'PARAMETERIZED'
       FRAME_-28976_FAMILY           = 'TWO-VECTOR'
       FRAME_-28976_PRI_AXIS         = 'Z'
       FRAME_-28976_PRI_VECTOR_DEF   = 'CONSTANT'
       FRAME_-28976_PRI_SPEC         = 'RECTANGULAR'
       FRAME_-28976_PRI_VECTOR       = ( 0, 0, 1 )
       FRAME_-28976_PRI_FRAME        = 'JUICE_JUPITER_MAG_SIII'       
       FRAME_-28976_SEC_AXIS         = 'X'
       FRAME_-28976_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
       FRAME_-28976_SEC_OBSERVER     = 'JUPITER'
       FRAME_-28976_SEC_TARGET       = 'SUN'
       FRAME_-28976_SEC_ABCORR       = 'NONE'

  \begintext


Jupiter Solar Wind frame (JUICE_JSW)
------------------------------------------------------------------------

   Definition:
   -----------  
   The Jupiter Solar Wind frame is defined is defined, based on the
   definition in [12], section 3.3.7) for the Geocentric Solar Wind 
   frame, as follows:
   
      -  +Z is perpendicular to the orbital plane of Jupiter and is 
         positive towards North;
         
      -  +X axis is the component of the solar wind direction vector
         that is orthogonal to the +Z-axis. This axis lies in the
         orbital plane of Jupiter and is positive in the direction 
         opposite to the solar wind;
         
      -  +Y axis completes the right-handed system;
         
      -  the origin of this frame is the centre of mass of
         Jupiter.

   
   Uses and applications:
   ----------------------
   The Jupiter Solar Wind frame is useful in analyzing the impact
   of solar wind on hemispheric events (from [12], section 3.3.7).
   
   
   Remarks:
   --------
   A critical issue to consider in the definition of this 
   frame is aberration – the bow shock of Jupiter is rotated in 
   the plane of the orbit of Jupiter by angle Vplanet/Vsolar_wind, 
   where Vplanet is the velocity of Jupiter in its orbit and 
   Vsolar_wind is the solar wind speed. The sense of rotation is 
   such the bow shock lags its un-rotated location on the anti-sunward 
   side of Jupiter. This is the +Y direction for Jupiter solar 
   orbital frame (JUPITER_SUN_ORB). Strictly speaking Vplanet should 
   be the component perpendicular to the solar wind. 
   
   Th solar wind streams off of the Sun in all directions at speeds 
   of about 400 km/s. Nevertheless, the solar wind is not uniform. 
   Although it is always directed away from the Sun, it changes speed 
   and carries with it magnetic clouds, interacting regions where high
   speed wind catches up with slow speed wind, and composition 
   variations. The solar wind speed can range from high (800 km/s) over
   coronal holes to low velocities (300 km/s) over streamers 
   (see [13]).
   
   For Jupiter under normal solar wind conditions (v ~400 km/s) the 
   angle of aberration ranges from 1.9652 degrees for the maximum 
   orbital velocity to 1.7819 degrees for the minimum orbital 
   velocity of Jupiter.
   
   By convention, the average velocity for Jupiter in its orbit and for
   the solar wind at Jupiter will be used in the definition of this
   frame, which is highly significant for determining  magnetopause and
   bow shock locations. This value is determined to be 1.8707 degrees. 
   
   Since the minimum aberration is of approximately 0.89095 degrees
   (maximum solar wind and minimum Jupiter orbital velocities) and 
   the maximum aberration is of approximately 2.6203 degrees (minimum 
   solar wind and maximum Jupiter orbital velocities), the error in 
   the definition of the +X-axis direction is between -0.97975 to 
   0.7496 degrees.
      
   This frame is defined relative to JUPITER_SUN_ORB, which is 
   dynamic. This aspect shall be taken into account if/when using this
   frame to define other dynamic frames. For further details, please
   refer to [1].
   
   
  \begindata

      FRAME_JUICE_JSW           =  -28977
      FRAME_-28977_NAME         = 'JUICE_JSW'
      FRAME_-28977_CLASS        =  4
      FRAME_-28977_CLASS_ID     =  -28977
      FRAME_-28977_CENTER       =  599
      TKFRAME_-28977_SPEC       = 'ANGLES'
      TKFRAME_-28977_RELATIVE   = 'JUPITER_SUN_ORB'
      TKFRAME_-28977_ANGLES     = (   -1.8707,  0.0,  0.0 )
      TKFRAME_-28977_AXES       = (         3,    2,    3 )
      TKFRAME_-28977_UNITS      = 'DEGREES'

  \begintext

   
Jupiter Solar Wind Magnetospheric frame (JUICE_JSWM)
------------------------------------------------------------------------

   Definition:
   -----------  
   The Jupiter solar wind magnetospheric frame is defined, based on the
   definition in [12], section 3.3.6) for the Geocentric Solar Wind 
   Magnetospheric frame, as follows:
   
      -  +X axis is aligned with the solar wind direction. This axis is
         positive in the direction opposite to the solar wind;
         
      -  +Z axis is the component of the magnetic axis (magnetic
         dipole) orthogonal to the +X axis, positive towards magnetic
         north of Jupiter;
         
      -  +Y axis completes the right-handed system. It is perpendicular
         to the magnetic axis (magnetic dipole);
         
      -  the origin of this frame is the centre of mass of the 
         Jupiter.

   All vectors are geometric: no corrections are used.
   
   
   Remarks:
   --------
   A critical issue to consider in the definition of this 
   frame is aberration – the bow shock of Jupiter is rotated in 
   the plane of the orbit of Jupiter by angle Vplanet/Vsolar_wind, 
   where Vplanet is the velocity of Jupiter in its orbit and 
   Vsolar_wind is the solar wind speed. The sense of rotation is 
   such the bow shock lags its un-rotated location on the anti-sunward 
   side of Jupiter. This is the +Y direction for Jupiter Solar 
   Orbital frame (JUPITER_SUN_ORB). Strictly speaking Vplanet should 
   be the component perpendicular to the solar wind. 
   
   Th solar wind streams off of the Sun in all directions at speeds 
   of about 400 km/s. Nevertheless, the solar wind is not uniform. 
   Although it is always directed away from the Sun, it changes speed 
   and carries with it magnetic clouds, interacting regions where high
   speed wind catches up with slow speed wind, and composition 
   variations. The solar wind speed can range from high (800 km/s) over
   coronal holes to low velocities (300 km/s) over streamers 
   (see [13]).
   
   For Jupiter under normal solar wind conditions (v ~400 km/s) the 
   angle of aberration ranges from 1.9652 degrees for the maximum 
   orbital velocity to 1.7819 degrees for the minimum orbital 
   velocity of Jupiter.
   
   By convention, the average velocity for Jupiter in its orbit and for
   the solar wind at Jupiter, will be used in the definition of this
   frame, which is highly significant for determining magnetopause and
   bow shock locations. This value is determined to be 1.8707 degrees. 
   
   Since the minimum aberration is of approximately 0.89095 degrees
   (maximum solar wind and minimum Jupiter orbital velocities) and 
   the maximum aberration is of approximately 2.6203 degrees (minimum 
   solar wind and maximum Jupiter orbital velocities), the error in 
   the definition of the +X-axis direction is between -0.97975 to 
   0.7496 degrees.
   
   This frame is defined relative to JUICE_JUPITER_BSM, which is 
   dynamic. This aspect shall be taken into account if/when using this
   frame to define other dynamic frames. For further details, please
   refer to [1].
   
  \begindata

      FRAME_JUICE_JSWM          =  -28978
      FRAME_-28978_NAME         = 'JUICE_JSWM'
      FRAME_-28978_CLASS        =  4
      FRAME_-28978_CLASS_ID     =  -28978
      FRAME_-28978_CENTER       =  599
      TKFRAME_-28978_SPEC       = 'ANGLES'
      TKFRAME_-28978_RELATIVE   = 'JUICE_JUPITER_BSM'
      TKFRAME_-28978_ANGLES     = (   -1.8707,  0.0,  0.0 )
      TKFRAME_-28978_AXES       = (         3,    2,    3 )
      TKFRAME_-28978_UNITS      = 'DEGREES'

  \begintext


Callisto Orbital frame transformed (JUICE_CALLISTO_PHI_ORB)
------------------------------------------------------------------------
   
   Definition:
   -----------
   The Callisto orbital frame is defined as follows:

      -  Y axis is the position of the Jupiter relative to 
         Callisto; it's the primary vector and points 
         from Callisto to Jupiter;
 
      -  -X axis is the component of the inertially reference
         velocity of Jupiter relative to Callisto orthogonal 
         to the +Y axis;
   
      -  +Z axis completes the right-handed system;
   
      -  the origin of this frame is the center of mass of
         Callisto.
   
   All vectors are geometric: no corrections are used.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Callisto vector in J2000 reference frame have 
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity' 
   vector. Therefore, the ephemeris data required to compute the 
   Callisto-Jupiter velocity vector in the J2000 reference frame 
   have to be loaded before using this frame.


   Remarks:
   --------
   This frame is a transformation of the CALLISTO_JUPITER_ORB frame.
   The transformation applied is x'=-y, y'=x, z'=z. 

   This frame is defined based on SPK data: different planetary 
   ephemerides for Callisto, Jupiter and the Jupiter Barycenter
   will lead to a different frame orientation at a given time.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, e.g. 
   JUICE_CALLISTO_PHI_ORB using de405 ephemerides.


  \begindata

      FRAME_JUICE_CALLISTO_PHI_ORB   =  -28980
      FRAME_-28980_NAME              = 'JUICE_CALLISTO_PHI_ORB' 
      FRAME_-28980_CLASS             =  5
      FRAME_-28980_CLASS_ID          =  -28980
      FRAME_-28980_CENTER            =  504
      FRAME_-28980_RELATIVE          = 'J2000'
      FRAME_-28980_DEF_STYLE         = 'PARAMETERIZED'
      FRAME_-28980_FAMILY            = 'TWO-VECTOR'
      FRAME_-28980_PRI_AXIS          = 'Y'
      FRAME_-28980_PRI_VECTOR_DEF    = 'OBSERVER_TARGET_POSITION'
      FRAME_-28980_PRI_OBSERVER      = 'CALLISTO'
      FRAME_-28980_PRI_TARGET        = 'JUPITER'
      FRAME_-28980_PRI_ABCORR        = 'NONE'
      FRAME_-28980_SEC_AXIS          = '-X'
      FRAME_-28980_SEC_VECTOR_DEF    = 'OBSERVER_TARGET_VELOCITY'
      FRAME_-28980_SEC_OBSERVER      = 'CALLISTO'
      FRAME_-28980_SEC_TARGET        = 'JUPITER'
      FRAME_-28980_SEC_ABCORR        = 'NONE'
      FRAME_-28980_SEC_FRAME         = 'J2000'
     
  \begintext
  

Europa Orbital frame transformed (JUICE_EUROPA_PHI_ORB)
------------------------------------------------------------------------
   
   Definition:
   -----------
   The Europa Orbital frame is defined as follows:

      -  +Y axis is the position of the Jupiter relative to 
         Europa; it's the primary vector and points 
         from Europa to Jupiter;
 
      -  -X axis is the component of the inertially referenced
         velocity of Jupiter relative to Europa orthogonal 
         to the +Y axis;

      -  +Z axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Europa.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Europa vector in J2000 reference frame have 
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity' 
   vector. Therefore, the ephemeris data required to compute the 
   Europa-Jupiter velocity vector in the J2000 reference frame 
   have to be loaded before using this frame.


   Remarks:
   --------
   This frame is a transformation of the EUROPA_JUPITER_ORB frame.
   The transformation applied is x'=-y, y'=x, z'=z. 

   This frame is defined based on SPK data: different planetary 
   ephemerides for Europa, Jupiter and the Jupiter Barycenter
   will lead to a different frame orientation at a given time.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, e.g. 
   JUICE_EUROPA_PHI_ORB using de405 ephemerides.


  \begindata

      FRAME_JUICE_EUROPA_PHI_ORB      =  -28981
      FRAME_-28981_NAME               = 'JUICE_EUROPA_PHI_ORB' 
      FRAME_-28981_CLASS              =  5
      FRAME_-28981_CLASS_ID           =  -28981
      FRAME_-28981_CENTER             =  502
      FRAME_-28981_RELATIVE           = 'J2000'
      FRAME_-28981_DEF_STYLE          = 'PARAMETERIZED'
      FRAME_-28981_FAMILY             = 'TWO-VECTOR'
      FRAME_-28981_PRI_AXIS           = 'Y'
      FRAME_-28981_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
      FRAME_-28981_PRI_OBSERVER       = 'EUROPA'
      FRAME_-28981_PRI_TARGET         = 'JUPITER'
      FRAME_-28981_PRI_ABCORR         = 'NONE'
      FRAME_-28981_SEC_AXIS           = '-X'
      FRAME_-28981_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
      FRAME_-28981_SEC_OBSERVER       = 'EUROPA'
      FRAME_-28981_SEC_TARGET         = 'JUPITER'
      FRAME_-28981_SEC_ABCORR         = 'NONE'
      FRAME_-28981_SEC_FRAME          = 'J2000'

  \begintext
  

Ganymede Orbital frame transformed (JUICE_GANYMEDE_PHI_ORB)
------------------------------------------------------------------------
   
   Definition:
   -----------
   The Ganymede Orbital frame is defined as follows:

      -  +Y axis is the position of the Jupiter relative to 
         Ganymede; it's the primary vector and points 
         from Ganymede to Jupiter;
 
      -  -X axis is the component of the inertially referenced
         velocity of Jupiter relative to Ganymede orthogonal 
         to the +Y axis;

      -  +Z axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Ganymede.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Ganymede vector in J2000 reference frame have 
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity' 
   vector. Therefore, the ephemeris data required to compute the 
   Ganymede-Jupiter velocity vector in the J2000 reference frame 
   have to be loaded before using this frame.


   Remarks:
   --------
   This frame is a transformation of the EUROPA_JUPITER_ORB frame.
   The transformation applied is x'=-y, y'=x, z'=z. 

   This frame is defined based on SPK data: different planetary 
   ephemerides for Ganymede, Jupiter and the Jupiter Barycenter
   will lead to a different frame orientation at a given time.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, e.g. 
   JUICE_GANYMEDE_PHI_ORB using de405 ephemerides.


  \begindata

      FRAME_JUICE_GANYMEDE_PHI_ORB    =  -28982
      FRAME_-28982_NAME               = 'JUICE_GANYMEDE_PHI_ORB'
      FRAME_-28982_CLASS              =  5
      FRAME_-28982_CLASS_ID           =  -28982
      FRAME_-28982_CENTER             =  503
      FRAME_-28982_RELATIVE           = 'J2000'
      FRAME_-28982_DEF_STYLE          = 'PARAMETERIZED'
      FRAME_-28982_FAMILY             = 'TWO-VECTOR'
      FRAME_-28982_PRI_AXIS           = 'Y'
      FRAME_-28982_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
      FRAME_-28982_PRI_OBSERVER       = 'GANYMEDE'
      FRAME_-28982_PRI_TARGET         = 'JUPITER'
      FRAME_-28982_PRI_ABCORR         = 'NONE'
      FRAME_-28982_SEC_AXIS           = '-X'
      FRAME_-28982_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
      FRAME_-28982_SEC_OBSERVER       = 'GANYMEDE'
      FRAME_-28982_SEC_TARGET         = 'JUPITER'
      FRAME_-28982_SEC_ABCORR         = 'NONE'
      FRAME_-28982_SEC_FRAME          = 'J2000'


  \begintext


JUICE Heliocentric Radial-Tangential-Normal (JUICE_HGRTN)
------------------------------------------------------------------------

   SPICE Frame Name and literature references
   ------------------------------------------
   Within the SPICE system, the JUICE Heliocentric
   Radial-Tangential-Normal frame is referred as JUICE_HGRTN. In
   literature, this frame is referred as HGRTN (from [4], [6], [7], and
   [8]).
   
      
   Definition:
   -----------
   The JUICE Heliocentric Radial-Tangential-Normal frame is defined as 
   follows (from [4]):

      -  the position of JUICE relative to the Sun is the 
         primary vector: +X axis points from the Sun to
         JUICE;

      -  the projection of the solar rotational axis perpendicular to
         the +X axis defines the frame's +Z axis;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the Sun's center of mass.

   All vectors are geometric: no corrections are used.
   
   
   Uses and applications (from [6]):
   ----------------------------------
   This frame is used to define the velocity and field direction of the 
   plasma environment that the spacecraft finds itself in.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different types
   of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' 
   vector, therefore, the ephemeris data required to compute the 
   JUICE-Sun vector in J2000 frame have to be loaded before 
   using this frame.

   The secondary vector is defined as a constant vector in the IAU_SUN
   frame, which is a PCK-based frame, therefore a PCK file containing
   the orientation constants for the Sun has to be loaded before using 
   this frame.


   Remarks:
   --------
   This frame is defined based on the IAU_SUN frame, whose evaluation is
   based on the data included in the loaded PCK file: different
   orientation constants for the Sun's spin axis will lead to different
   frames.

   Since the primary vector of this frame is defined as an
   'observer-target position' vector, the usage of different ephemerides 
   for the Sun, the Solar System Barycentre and JUICE may lead to 
   different orientations of this frame.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, i.e. 
   JUICE_HGRTN using IAU 2009 constants, de405 ephemerides family and 
   JUICE ephemerides version N.


  \begindata

      FRAME_JUICE_HGRTN             = -28990
      FRAME_-28990_NAME             = 'JUICE_HGRTN' 
      FRAME_-28990_CLASS            =  5
      FRAME_-28990_CLASS_ID         =  -28990
      FRAME_-28990_CENTER           =  10
      FRAME_-28990_RELATIVE         = 'J2000'
      FRAME_-28990_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-28990_FAMILY           = 'TWO-VECTOR'
      FRAME_-28990_PRI_AXIS         = 'X'
      FRAME_-28990_PRI_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-28990_PRI_OBSERVER     = 'SUN'
      FRAME_-28990_PRI_TARGET       = 'JUICE'
      FRAME_-28990_PRI_ABCORR       = 'NONE'
      FRAME_-28990_SEC_AXIS         = 'Z'
      FRAME_-28990_SEC_VECTOR_DEF   = 'CONSTANT'
      FRAME_-28990_SEC_FRAME        = 'IAU_SUN'
      FRAME_-28990_SEC_SPEC         = 'RECTANGULAR'
      FRAME_-28990_SEC_VECTOR       = ( 0, 0, 1 )

  \begintext


Sun JUICE Radial-Tangential-Normal (JUICE_SUN_RTN)
------------------------------------------------------------------------

   SPICE frame name and literature references
   ------------------------------------------
   Within the SPICE system, the Sun JUICE Radial-Tangential-Normal
   frame is referred as JUICE_SUN_RTN. In literature, this frame is
   referred as RTN (from [4]), Orbit RTN Coordinate System (from 
   [10]), or radial-transverse-normal (from [5]).
   
   
   Definition:
   -----------
   The Sun JUICE Radial-Tangential-Normal frame is defined as
   follows (from [10]):
   
      -  the position of Sun relative to JUICE is the 
         primary vector: +X axis points from the JUICE to
         Sun;

      -  the normal vector to the orbit plane defines the 
         frame's +Z axis;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is JUICE's center of mass.

   All vectors are geometric: no corrections are used.
   
   
   Uses and applications (from [10]):
   ----------------------------------
   This frame is used to project the perturbation forces imposed on the
   spacecraft so that the Gauss Lagrange equations can be established
   and the variations of the orbit elements expressed as function
   of the three components of the perturbation acceleration in the
   orbit radial (+X), tangential (+Y) and normal (+Z) axis.
   
   
   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different types
   of specifications for the primary and secondary vectors.
   
   The primary vector is defined as an 'observer-target position' 
   vector and the secondary vector is defined as an 'observer-target'
   velocity, therefore, the ephemeris data required to compute the 
   Sun-JUICE state in J2000 frame have to be loaded before 
   using this frame.
   
   
   Remarks:
   --------
   Since the primary and secondary vectors of this frame are defined
   based on the Sun-JUICE state vector, the usage of 
   different ephemerides for the Sun, Solar System Barycentre and 
   for JUICE conduce also to different frame orientation of this frame
   at a given time.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, i.e. 
   JUICE_SUN_RTN using de405 ephemerides family and JUICE ephemerides
   version N.
   
   
   \begindata

      FRAME_JUICE_SUN_RTN           = -28991
      FRAME_-28991_NAME             = 'JUICE_SUN_RTN'
      FRAME_-28991_CLASS            =  5
      FRAME_-28991_CLASS_ID         = -28991
      FRAME_-28991_CENTER           = -28
      FRAME_-28991_RELATIVE         = 'J2000'
      FRAME_-28991_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-28991_FAMILY           = 'TWO-VECTOR'
      FRAME_-28991_PRI_AXIS         = 'X'
      FRAME_-28991_PRI_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-28991_PRI_OBSERVER     = 'JUICE'
      FRAME_-28991_PRI_TARGET       = 'SUN'
      FRAME_-28991_PRI_ABCORR       = 'NONE'
      FRAME_-28991_SEC_AXIS         = 'Y'
      FRAME_-28991_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_VELOCITY'
      FRAME_-28991_SEC_OBSERVER     = 'JUICE'
      FRAME_-28991_SEC_TARGET       = 'SUN'
      FRAME_-28991_SEC_ABCORR       = 'NONE'
      FRAME_-28991_SEC_FRAME        = 'J2000'

   \begintext     


Jupiter JUICE Dipole Meridian (JUICE_JUPITER_DM)
------------------------------------------------------------------------

   SPICE Frame Name and literature references
   ------------------------------------------
   Within the SPICE system, the Jupiter JUICE Dipole Meridian Frame 
   frame is referred as JUICE_JUPITER_DM. In literature, this frame 
   is referred as Jupiter DM (from [3] and [11]).


   Definition:
   -----------
   The Jupiter JUICE Dipole Meridian frame is defined, based on the definition
   provided in [11] for the Earth's dipole meridian frame, as follows:
   
      -  +Z axis is aligned with Jupiter's north magnetic pole;

      -  the Y-axis is chosen to be perpendicular to a radius vector to 
         the point of observation. The positive Y direction is chosen 
         to be eastwards, so that the X-axis is directed outwards from 
         Jupiter to JUICE.

      -  the X-axis completes the right-handed system;

      -  the origin of this frame is Jupiter's center of mass.
     
   All vectors are geometric: no corrections are used.
   
   
   Uses and applications (from [11]):
   ----------------------------------
   It is used to order data controlled by the dipole magnetic field
   where the influence of the solar wind interaction with the
   magnetosphere is weak. It has been used extensively to describe the
   distortions of the magnetospheric field in terms of the two angles
   declination and inclination which can be easily derived from
   measurements in this frame.
   

   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different types
   of specifications for the primary and secondary vectors.
   
   The primary vector is defined as a constant vector in the 
   JUICE_JUPITER_MAG body-fixed frame, which provides the orientation
   for the Jupiter North magnetic dipole and the Jupiter magnetic
   frame with respect to the IAU_JUPITER frame.   
   
   The secondary vector is defined as an 'observer-target position' 
   vector, therefore, the ephemeris data required to compute the 
   JUICE-Jupiter vector in J2000 frame have to be loaded before 
   using this frame.


   Remarks:
   --------
   The primary vector of this frame is defined based on the 
   JUICE_JUPITER_MAG frame. Different implementations of this frame, i.e. 
   changes to the existing frame definition will lead to different
   orientations of the North geomagnetic centered dipole, but only when
   these data lead to different right ascension and declination
   coordinates of the North magnetic centered dipole vector of Jupiter.
        
   Since the secondary vector of this frame is defined as an 
   'observer-target position' vector, the usage of different 
   ephemerides for JUICE, Jupiter and Jupiter barycentre may lead to 
   different orientations of this frame, but only when these data lead to 
   different projections of the Jupiter-JUICE vector on the magnetic 
   equator of Jupiter.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, i.e. 
   JUICE_JUPITER_DM using IAU 2009 constants and de403 ephemerides.
   

  \begindata

      FRAME_JUICE_JUPITER_DM        = -28992
      FRAME_-28992_NAME             = 'JUICE_JUPITER_DM'
      FRAME_-28992_CLASS            =  5
      FRAME_-28992_CLASS_ID         = -28992
      FRAME_-28992_CENTER           =  599
      FRAME_-28992_RELATIVE         = 'J2000'
      FRAME_-28992_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-28992_FAMILY           = 'TWO-VECTOR'
      FRAME_-28992_PRI_AXIS         = 'Z'
      FRAME_-28992_PRI_VECTOR_DEF   = 'CONSTANT'
      FRAME_-28992_PRI_FRAME        = 'JUICE_JUPITER_MAG'
      FRAME_-28992_PRI_SPEC         = 'RECTANGULAR'
      FRAME_-28992_PRI_VECTOR       = ( 0, 0, 1 )
      FRAME_-28992_SEC_AXIS         = 'X'
      FRAME_-28992_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-28992_SEC_OBSERVER     = 'JUPITER'
      FRAME_-28992_SEC_TARGET       = 'JUICE'
      FRAME_-28992_SEC_ABCORR       = 'NONE'

  \begintext


JUPITER Generic Frame Definitions
-------------------------------------------------------------------------------

   This section contains the definition of the Jupiter System generic frames.
  

Jupiter Mean Equator of date frame (JUPITER_MEQUD)
------------------------------------------------------------------------

   Definition:
   -----------   
   The Jupiter Mean Equator of date frame is defined as follows:

      -  X-Y plane is defined by Jupiter equator of date, and
         the +Z axis is parallel to Jupiter rotation axis 
         of date, pointing toward the North side of the 
         invariable plane;

      -  +X axis is the component of the ascending node of Jupiter
         equator of date on the Earth Mean Equator of J2000 
         orthogonal to the +Z axis;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of 
         Jupiter.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using constant vectors 
   as the specification method.

   The primary vector is defined as a constant vector in the 
   IAU_JUPITER frame, which is a PCK-based frame. Therefore, a PCK
   file containing the orientation constants for Jupiter has to be 
   loaded before using this frame.
   
   The secondary vector is defined in the J2000 frame and therefore
   it does not require to load any additional data.


   Remarks:
   --------
   This frame is defined based on the IAU_JUPITER frame, whose 
   evaluation is based on the data included in the loaded PCK file: 
   different orientation constants for the spin axis of Jupiter will
   lead to a different frame orientation at a given time.

   This frame is provided as the ``most generic'' Jupiter mean 
   equator of date frame since the user has the possibility of loading 
   different Jupiter orientation constants that would help to 
   define different implementations of this frame.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, e.g. 
   JUPITER_MEQUD using IAU 2009 constants.


  \begindata

      FRAME_JUPITER_MEQUD             =  500599000
      FRAME_500599000_NAME            = 'JUPITER_MEQUD' 
      FRAME_500599000_CLASS           =  5
      FRAME_500599000_CLASS_ID        =  500599000
      FRAME_500599000_CENTER          =  599
      FRAME_500599000_RELATIVE        = 'J2000'
      FRAME_500599000_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500599000_FAMILY          = 'TWO-VECTOR'
      FRAME_500599000_PRI_AXIS        = 'Z'
      FRAME_500599000_PRI_VECTOR_DEF  = 'CONSTANT'
      FRAME_500599000_PRI_FRAME       = 'IAU_JUPITER' 
      FRAME_500599000_PRI_SPEC        = 'RECTANGULAR'
      FRAME_500599000_PRI_VECTOR      = ( 0, 0, 1 )
      FRAME_500599000_SEC_AXIS        = 'Y'
      FRAME_500599000_SEC_VECTOR_DEF  = 'CONSTANT'                
      FRAME_500599000_SEC_FRAME       = 'J2000'
      FRAME_500599000_SEC_SPEC        = 'RECTANGULAR'
      FRAME_500599000_SEC_VECTOR      = ( 0, 0, 1 )

  \begintext


Jupiter Solar Equatorial frame (JUPITER_SUN_EQU)
------------------------------------------------------------------------

   Definition:
   -----------   
   The Jupiter Solar Equatorial frame is defined, based on the
   definition provided in [15] for the solar equatorial
   frame, as follows:

      -  +X axis is the position of the Sun relative to Jupiter; 
         it's the primary vector and points from Jupiter to the 
         Sun;

      -  +Z axis is the component of the Sun's north pole of date 
         orthogonal to the +X axis;

      -  +Y axis completes the right-handed reference frame;

      -  the origin of this frame is center of mass of Jupiter.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.
   
   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Sun vector in J2000 reference frame have to be loaded 
   before using this frame.
   
   The primary vector is defined as a constant vector in the IAU_SUN
   frame, which is a PCK-based frame. Therefore a PCK file containing
   the orientation constants for Sun has to be loaded before using 
   this frame.


   Remarks:
   --------
   This frame is defined based on SPK data: different planetary
   ephemerides for Jupiter, Jupiter Barycenter, the Sun and the Solar 
   System Barycenter will lead to different frame orientation at a 
   given time.
   
   This frame is also defined based on the IAU_SUN frame, whose 
   evaluation is based on the data included in the loaded PCK file: 
   different orientation constants for the Sun's spin axis will lead 
   to different frame orientation at a given time. 
   
   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, i.e. 
   JUPITER_SUN_EQU using IAU 2009 constants and de403 ephemerides.


  \begindata

      FRAME_JUPITER_SUN_EQU           =  500599001
      FRAME_500599001_NAME            = 'JUPITER_SUN_EQU' 
      FRAME_500599001_CLASS           =  5
      FRAME_500599001_CLASS_ID        =  500599001
      FRAME_500599001_CENTER          =  599
      FRAME_500599001_RELATIVE        = 'J2000'
      FRAME_500599001_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500599001_FAMILY          = 'TWO-VECTOR'
      FRAME_500599001_PRI_AXIS        = 'X'
      FRAME_500599001_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_500599001_PRI_OBSERVER    = 'JUPITER'
      FRAME_500599001_PRI_TARGET      = 'SUN'
      FRAME_500599001_PRI_ABCORR      = 'NONE'
      FRAME_500599001_SEC_AXIS        = 'Z'
      FRAME_500599001_SEC_VECTOR_DEF  = 'CONSTANT'
      FRAME_500599001_SEC_FRAME       = 'IAU_SUN' 
      FRAME_500599001_SEC_SPEC        = 'RECTANGULAR'
      FRAME_500599001_SEC_VECTOR      = ( 0, 0, 1 )


  \begintext


Jupiter Solar Orbital frame (JUPITER_SUN_ORB)
------------------------------------------------------------------------
   
   SPICE frame name and literature references
   ------------------------------------------
   Within the SPICE system, the Jupiter Solar Orbital is referred as
   JUPITER_SUN_ORB. In literature (see [11]) and in the documentation
   provided by the JUICE MAG team (see [14]), this frame is referred
   as JSO.

   
   Definition:
   -----------
   The Jupiter Solar Orbital frame is defined as follows (from [14]):

      -  +X axis is the position of the Sun relative to 
         Jupiter; it's the primary vector and points 
         from Jupiter to Sun;
 
      -  +Y axis is the component of the inertially referenced
         velocity of Sun relative to Jupiter orthogonal 
         to the +X axis;

      -  +Z axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Jupiter.

   All vectors are geometric: no corrections are used.
   
   
   Uses and applications:
   ----------------------
   The Jupiter Solar Orbital frame is used primarily for orbit analysis,
   the analysis of solar wind data and interplanetary magnetic field data,
   and other data where the solar wind flow direction (-x) is the dominant
   process.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Sun vector in J2000 reference frame have 
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Sun velocity vector in the J2000 reference frame 
   have to be loaded before using this frame.


   Remarks:
   --------
   This frame is defined based on SPK data: different planetary 
   ephemerides for Jupiter, Sun and the Sun Barycenter
   will lead to a different frame orientation at a given time.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, e.g. 
   JUPITER_SUN_ORB using de405 ephemerides.


  \begindata

      FRAME_JUPITER_SUN_ORB           =  500599002
      FRAME_500599002_NAME            = 'JUPITER_SUN_ORB' 
      FRAME_500599002_CLASS           =  5
      FRAME_500599002_CLASS_ID        =  500599002
      FRAME_500599002_CENTER          =  599
      FRAME_500599002_RELATIVE        = 'J2000'
      FRAME_500599002_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500599002_FAMILY          = 'TWO-VECTOR'
      FRAME_500599002_PRI_AXIS        = 'X'
      FRAME_500599002_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_500599002_PRI_OBSERVER    = 'JUPITER'
      FRAME_500599002_PRI_TARGET      = 'SUN'
      FRAME_500599002_PRI_ABCORR      = 'NONE'
      FRAME_500599002_SEC_AXIS        = 'Y'
      FRAME_500599002_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
      FRAME_500599002_SEC_OBSERVER    = 'JUPITER'
      FRAME_500599002_SEC_TARGET      = 'SUN'
      FRAME_500599002_SEC_ABCORR      = 'NONE'
      FRAME_500599002_SEC_FRAME       = 'J2000'

  \begintext


Jupiter-centric Callisto-following frame (JUPITER_CALLISTO_BCSF)
------------------------------------------------------------------------

   Definition:
   -----------
   The Jupiter-centric Callisto-following frame is defined as 
   follows:
            
	  -  X-Y plane is defined by the Jupiter's orbital plane 
	     of date: the +Z axis, primary vector, is the normal 
	     vector to this plane, always pointing toward the North 
	     side of the invariable plane;

	  -  +X axis is the component of the Jupiter-Callisto 
	     vector that is orthogonal to the +Z axis;

	  -  +Y axis completes the right-handed system;

	  -  the origin of this frame is the center of mass of the
	     Jupiter
	  
   All vectors are geometric: no corrections are used.


   Required Data:
   --------------

   The frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as a constant vector in the 
   JUPITER_SUN_ORB (Jupiter orbital) frame, which is a dynamic 
   frame. Therefore, the data required to evaluate that base frame 
   at the requested epoch needs to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Callisto vector in J2000 frame have to be loaded
   before using this frame.
   
   
   Remarks:
   --------
   SPICE imposes a constraint in the definition of dynamic frames 
   (see [1]):

      When the definition of a parameterized dynamic frame F1 refers to
      a second frame F2 the referenced frame F2 may be dynamic, but F2
      must not make reference to any dynamic frame.

   Therefore, no other dynamic frame should make reference to this
   frame.
   
   This frame is defined based on SPK data: different planetary 
   ephemerides for Callisto, Jupiter, Jupiter Barycenter, the Sun 
   and the Solar System will lead to a different frame orientation 
   at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g. 
   JUPITER_CALLISTO_BCSF using de405 ephemerides. 

   
  \begindata

      FRAME_JUPITER_CALLISTO_BCSF     =  500599006
      FRAME_500599006_NAME            = 'JUPITER_CALLISTO_BCSF'
      FRAME_500599006_CLASS           =  5
      FRAME_500599006_CLASS_ID        =  500599006
      FRAME_500599006_CENTER          =  599
      FRAME_500599006_RELATIVE        = 'J2000'
      FRAME_500599006_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500599006_FAMILY          = 'TWO-VECTOR'
      FRAME_500599006_PRI_AXIS        = 'Z'
      FRAME_500599006_PRI_VECTOR_DEF  = 'CONSTANT'
      FRAME_500599006_PRI_FRAME       = 'JUPITER_SUN_ORB'
      FRAME_500599006_PRI_SPEC        = 'RECTANGULAR'
      FRAME_500599006_PRI_VECTOR      = ( 0, 0, 1 )
      FRAME_500599006_SEC_AXIS        = 'X'
      FRAME_500599006_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_500599006_SEC_OBSERVER    = 'JUPITER'
      FRAME_500599006_SEC_TARGET      = 'CALLISTO'  
      FRAME_500599006_SEC_ABCORR      = 'NONE'

  \begintext


Jupiter-centric Europa-following frame (JUPITER_EUROPA_BCSF)
------------------------------------------------------------------------

   Definition:
   -----------
   The Jupiter-centric Europa-following frame is defined as 
   follows:
            
	  -  X-Y plane is defined by the Jupiter's orbital plane 
	     of date: the +Z axis, primary vector, is the normal 
	     vector to this plane, always pointing toward the North 
	     side of the invariable plane;

	  -  +X axis is the component of the Jupiter-Europa 
	     vector that is orthogonal to the +Z axis;

	  -  +Y axis completes the right-handed system;

	  -  the origin of this frame is the center of mass of the
	     Jupiter
	  
   All vectors are geometric: no corrections are used.


   Required Data:
   --------------

   The frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as a constant vector in the 
   JUPITER_SUN_ORB (Jupiter orbital) frame, which is a dynamic 
   frame. Therefore, the data required to evaluate that base frame 
   at the requested epoch needs to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Europa vector in J2000 frame have to be loaded
   before using this frame.
   
   
   Remarks:
   --------
   SPICE imposes a constraint in the definition of dynamic frames 
   (see [1]):

      When the definition of a parameterized dynamic frame F1 refers to
      a second frame F2 the referenced frame F2 may be dynamic, but F2
      must not make reference to any dynamic frame.

   Therefore, no other dynamic frame should make reference to this
   frame.
   
   This frame is defined based on SPK data: different planetary 
   ephemerides for Europa, Jupiter, Jupiter Barycenter, the Sun 
   and the Solar System will lead to a different frame orientation 
   at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g. 
   JUPITER_EUROPA_BCSF using de405 ephemerides. 

   
  \begindata

      FRAME_JUPITER_EUROPA_BCSF       =  500599004
      FRAME_500599004_NAME            = 'JUPITER_EUROPA_BCSF'
      FRAME_500599004_CLASS           =  5
      FRAME_500599004_CLASS_ID        =  500599004
      FRAME_500599004_CENTER          =  599
      FRAME_500599004_RELATIVE        = 'J2000'
      FRAME_500599004_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500599004_FAMILY          = 'TWO-VECTOR'
      FRAME_500599004_PRI_AXIS        = 'Z'
      FRAME_500599004_PRI_VECTOR_DEF  = 'CONSTANT'
      FRAME_500599004_PRI_FRAME       = 'JUPITER_SUN_ORB'
      FRAME_500599004_PRI_SPEC        = 'RECTANGULAR'
      FRAME_500599004_PRI_VECTOR      = ( 0, 0, 1 )
      FRAME_500599004_SEC_AXIS        = 'X'
      FRAME_500599004_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_500599004_SEC_OBSERVER    = 'JUPITER'
      FRAME_500599004_SEC_TARGET      = 'EUROPA'  
      FRAME_500599004_SEC_ABCORR      = 'NONE'

  \begintext


Jupiter-centric Ganymede-following frame (JUPITER_GANYMEDE_BCSF)
------------------------------------------------------------------------

   Definition:
   -----------
   The Jupiter-centric Ganymede-following frame is defined as 
   follows:
            
	  -  X-Y plane is defined by the Jupiter's orbital plane 
	     of date: the +Z axis, primary vector, is the normal 
	     vector to this plane, always pointing toward the North 
	     side of the invariable plane;

	  -  +X axis is the component of the Jupiter-Ganymede 
	     vector that is orthogonal to the +Z axis;

	  -  +Y axis completes the right-handed system;

	  -  the origin of this frame is the center of mass of the
	     Jupiter
	  
   All vectors are geometric: no corrections are used.


   Required Data:
   --------------

   The frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as a constant vector in the 
   JUPITER_SUN_ORB (Jupiter orbital) frame, which is a dynamic 
   frame. Therefore, the data required to evaluate that base frame 
   at the requested epoch needs to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Jupiter-Ganymede vector in J2000 frame have to be loaded
   before using this frame.
   
   
   Remarks:
   --------
   SPICE imposes a constraint in the definition of dynamic frames 
   (see [1]):

      When the definition of a parameterized dynamic frame F1 refers to
      a second frame F2 the referenced frame F2 may be dynamic, but F2
      must not make reference to any dynamic frame.

   Therefore, no other dynamic frame should make reference to this
   frame.
   
   This frame is defined based on SPK data: different planetary 
   ephemerides for Ganymede, Jupiter, Jupiter Barycenter, the Sun 
   and the Solar System will lead to a different frame orientation 
   at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g. 
   JUPITER_GANYMEDE_BCSF using de405 ephemerides. 

   
  \begindata

      FRAME_JUPITER_GANYMEDE_BCSF     =  500599005
      FRAME_500599005_NAME            = 'JUPITER_GANYMEDE_BCSF'
      FRAME_500599005_CLASS           =  5
      FRAME_500599005_CLASS_ID        =  500599005
      FRAME_500599005_CENTER          =  599
      FRAME_500599005_RELATIVE        = 'J2000'
      FRAME_500599005_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500599005_FAMILY          = 'TWO-VECTOR'
      FRAME_500599005_PRI_AXIS        = 'Z'
      FRAME_500599005_PRI_VECTOR_DEF  = 'CONSTANT'
      FRAME_500599005_PRI_FRAME       = 'JUPITER_SUN_ORB'
      FRAME_500599005_PRI_SPEC        = 'RECTANGULAR'
      FRAME_500599005_PRI_VECTOR      = ( 0, 0, 1 )
      FRAME_500599005_SEC_AXIS        = 'X'
      FRAME_500599005_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_500599005_SEC_OBSERVER    = 'JUPITER'
      FRAME_500599005_SEC_TARGET      = 'GANYMEDE'  
      FRAME_500599005_SEC_ABCORR      = 'NONE'

  \begintext


Callisto Orbital frame (CALLISTO_JUPITER_ORB)
------------------------------------------------------------------------
   
   Definition:
   -----------
   The Callisto orbital frame is defined as follows:

      -  +X axis is the position of the Jupiter relative to 
         Callisto; it's the primary vector and points 
         from Callisto to Jupiter;
 
      -  +Y axis is the component of the inertially referenced
         velocity of Jupiter relative to Callisto orthogonal 
         to the +X axis;

      -  +Z axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Callisto.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Callisto-Jupiter vector in J2000 reference frame have 
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity' 
   vector. Therefore, the ephemeris data required to compute the 
   Callisto-Jupiter velocity vector in the J2000 reference frame 
   have to be loaded before using this frame.


   Remarks:
   --------
   This frame is defined based on SPK data: different planetary 
   ephemerides for Callisto, Jupiter and the Jupiter Barycenter
   will lead to a different frame orientation at a given time.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, e.g. 
   CALLISTO_JUPITER_ORB using de405 ephemerides.


  \begindata

      FRAME_CALLISTO_JUPITER_ORB      =  500504000
      FRAME_500504000_NAME            = 'CALLISTO_JUPITER_ORB' 
      FRAME_500504000_CLASS           =  5
      FRAME_500504000_CLASS_ID        =  500504000
      FRAME_500504000_CENTER          =  504
      FRAME_500504000_RELATIVE        = 'J2000'
      FRAME_500504000_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500504000_FAMILY          = 'TWO-VECTOR'
      FRAME_500504000_PRI_AXIS        = 'X'
      FRAME_500504000_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_500504000_PRI_OBSERVER    = 'CALLISTO'
      FRAME_500504000_PRI_TARGET      = 'JUPITER'
      FRAME_500504000_PRI_ABCORR      = 'NONE'
      FRAME_500504000_SEC_AXIS        = 'Y'
      FRAME_500504000_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
      FRAME_500504000_SEC_OBSERVER    = 'CALLISTO'
      FRAME_500504000_SEC_TARGET      = 'JUPITER'
      FRAME_500504000_SEC_ABCORR      = 'NONE'
      FRAME_500504000_SEC_FRAME       = 'J2000'

  \begintext
  

Europa Orbital frame (EUROPA_JUPITER_ORB)
------------------------------------------------------------------------
   
   Definition:
   -----------
   The Europa Orbital frame is defined as follows:

      -  +X axis is the position of the Jupiter relative to 
         Europa; it's the primary vector and points 
         from Europa to Jupiter;
 
      -  +Y axis is the component of the inertially referenced
         velocity of Jupiter relative to Europa orthogonal 
         to the +X axis;

      -  +Z axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Europa.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Europa-Jupiter vector in J2000 reference frame have 
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity' 
   vector. Therefore, the ephemeris data required to compute the 
   Europa-Jupiter velocity vector in the J2000 reference frame 
   have to be loaded before using this frame.


   Remarks:
   --------
   This frame is defined based on SPK data: different planetary 
   ephemerides for Europa, Jupiter and the Jupiter Barycenter
   will lead to a different frame orientation at a given time.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, e.g. 
   EUROPA_JUPITER_ORB using de405 ephemerides.


  \begindata

      FRAME_EUROPA_JUPITER_ORB        =  500502000
      FRAME_500502000_NAME            = 'EUROPA_JUPITER_ORB' 
      FRAME_500502000_CLASS           =  5
      FRAME_500502000_CLASS_ID        =  500502000
      FRAME_500502000_CENTER          =  502
      FRAME_500502000_RELATIVE        = 'J2000'
      FRAME_500502000_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500502000_FAMILY          = 'TWO-VECTOR'
      FRAME_500502000_PRI_AXIS        = 'X'
      FRAME_500502000_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_500502000_PRI_OBSERVER    = 'EUROPA'
      FRAME_500502000_PRI_TARGET      = 'JUPITER'
      FRAME_500502000_PRI_ABCORR      = 'NONE'
      FRAME_500502000_SEC_AXIS        = 'Y'
      FRAME_500502000_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
      FRAME_500502000_SEC_OBSERVER    = 'EUROPA'
      FRAME_500502000_SEC_TARGET      = 'JUPITER'
      FRAME_500502000_SEC_ABCORR      = 'NONE'
      FRAME_500502000_SEC_FRAME       = 'J2000'

  \begintext


Ganymede Orbital frame (GANYMEDE_JUPITER_ORB)
------------------------------------------------------------------------
   
   Definition:
   -----------
   The Ganymede Orbital frame is defined as follows:

      -  +X axis is the position of the Jupiter relative to 
         Ganymede; it's the primary vector and points 
         from Ganymede to Jupiter;
 
      -  +Y axis is the component of the inertially referenced
         velocity of Jupiter relative to Ganymede orthogonal 
         to the +X axis;

      -  +Z axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Ganymede.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------
   This frame is defined as a two-vector frame using two different 
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position' 
   vector. Therefore, the ephemeris data required to compute the 
   Ganymede-Jupiter vector in J2000 reference frame have 
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity' 
   vector. Therefore, the ephemeris data required to compute the 
   Ganymede-Jupiter velocity vector in the J2000 reference frame 
   have to be loaded before using this frame.


   Remarks:
   --------
   This frame is defined based on SPK data: different planetary 
   ephemerides for Ganymede, Jupiter and the Jupiter Barycenter
   will lead to a different frame orientation at a given time.
   
   It is strongly recommended to indicate what data have been used 
   in the evaluation of this frame when referring to it, e.g. 
   GANYMEDE_JUPITER_ORB using de405 ephemerides.


  \begindata

      FRAME_GANYMEDE_JUPITER_ORB      =  500503000
      FRAME_500503000_NAME            = 'GANYMEDE_JUPITER_ORB' 
      FRAME_500503000_CLASS           =  5
      FRAME_500503000_CLASS_ID        =  500503000
      FRAME_500503000_CENTER          =  503
      FRAME_500503000_RELATIVE        = 'J2000'
      FRAME_500503000_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_500503000_FAMILY          = 'TWO-VECTOR'
      FRAME_500503000_PRI_AXIS        = 'X'
      FRAME_500503000_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_500503000_PRI_OBSERVER    = 'GANYMEDE'
      FRAME_500503000_PRI_TARGET      = 'JUPITER'
      FRAME_500503000_PRI_ABCORR      = 'NONE'
      FRAME_500503000_SEC_AXIS        = 'Y'
      FRAME_500503000_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
      FRAME_500503000_SEC_OBSERVER    = 'GANYMEDE'
      FRAME_500503000_SEC_TARGET      = 'JUPITER'
      FRAME_500503000_SEC_ABCORR      = 'NONE'
      FRAME_500503000_SEC_FRAME       = 'J2000'

  \begintext

