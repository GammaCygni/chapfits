;---------------------------------------------------------------------------------------------------------
;Katy Fallows
;Boston University CSP
;
;GOAL: This procedure fits a double Chapman function to all input profiles and compiles the bestfits and 
;      fits arrays with all of the fit parameters and other information about each profile.
;      FOR MGS profiles, used to create save file mgschapfits.sav, which includes the array bestfits 
;      and bestfits_label, as well as good and bad fit definitions
;      Also created mgschapfits_gb.sav, which includes good/bad fit definitions

;5/24/2012: Creating fitsfn.pro from fits.pro - this version will work as a function, to work with any 
;           input electron densities. fits.pro was hard-wired to do only the MGS profiles.
;1/21/2013: copied fitsfn_roundtwo_cm3.pro    
;           making cm^-3 the final units; using two round fit; adding Gaussian bad fit definitions.
;           Hoping this will be the last time that I mess with this!       
;1/30/2013: Profiles with bad fit to N1 on the first round stay flagged as bad.  If not, can end up with 
;           fit that is marked as good but is very much not good
;3/20/2013: Added a different set of conditions for badfits = 1              
;---------------------------------------------------------------------------------------------------------

PRO fitsfn_final, eldenarr, altsarr, $                                                            ;inputs
        PH = ph, PL = pl, R2PL = r2pl, AIN = Ain, ZEN = zenarray, SAVEFILE = savefile, BADFITS = badfits, psind = psind, $  ;opt. inputs
        SHOWPLOTS = showplots, ISTART = istart, REV = rev, TOL = tol, ROUNDTWO = roundtwo, $      ;switches
        bestfits, bestfits_label, bestfits_comments, good, status                                 ;outputs
            
;-------------------------------------------INPUTS--------------------------------------------------------
;eldenarr  - number density array, n(x,y), where x = profiles, y = altitudes, in cm^-3
;altsarr   - altitude array, same dimensions as n, in km
;ph        - high altitude cut-off for fit - percentage of peak number density
;pl        - low altitude cut-off for fit, first round - percentage of peak number density 
;r2pl      - low altitude cut-off for the fit, second round - percentage of M1 peak density
;Ain       - initial guess of fit parameters; altitudes, for example, change enough that this matters
;zen       - optional array of solar zenith angles corresponding to profiles, just for display on plots
;savefile  - name of save file to be created, nothing saved if not specified
;SHOWPLOTS - if SHOWPLOTS is set, then each profile and fit are displayed in window one at a time
;psind     - index of profile plotted as .ps as an example
;ISTART    - index of first profile to start fitting - helpful if looking at plots using SHOWPLOTS
;REV       - NOT RECOMMENDED; reverse the arrays- if on, arrays descend with alt, matching original fit 
;TOL       - tolerance level for CURVEFIT; default = 1E-3
;ROUNDTWO  - RCOMMENDED; set to use two rounds of fitting, the second time using r2pl as the lower cutoff 
;BADFITS   - 0 - Fit Gaussian to distribution of each parameter, flag as bad those > N sigma from the mean
;          - 1 - Cut off points above or below hard cutoff limits (default)
;          - 2 - Alternate set of hard cutoff limits
;-------------------------------------------OUTPUTS-------------------------------------------------------
;bestfits          - array with bestfit parameters and some other useful info
;bestfits_label    - short title for each column in bestfits
;bestfits_comments - longer explanation of each column in bestfits
;status            - array of which bad fit flags for criterion and each profile
;---------------------------------------------------------------------------------------------------------

;testi = 3881    ;index of profile to follow for testing, set to -1 to skip testing lines
testi = -1

RESTORE, 'C:\Users\kfallows\IDLWorkspaces\General\hexcolors.sav'

jdateswitch = 'no'

nprof = N_ELEMENTS(eldenarr[*,0])
nparam = 6
bestfits=DBLARR(18,nprof)

bestfits_label    = [' 0 Profile_Number',' 1 Nmax',' 2 N_at_Zmax',' 3 N2',' 4 Z2',' 5 H2',' 6 N1',' 7 Z1',' 8 H1',' 9 Chi_Squared',$
                     '10 Status','11 SZA','12 n0','13 n1','14 z0','15 z1', '16 proRMS', '17 perRMS']
bestfits_comments = [' 0 Profile Number - a running index; for MGS Data, this corresponds to the index in all_eds_profiles', $
                     ' 1 Nmax - Maximum number density of the profile', $
                     ' 2 N at Zmax - Number density at the fitted height of M2', $
                     ' 3 N2 - Fit parameter, peak number density of M2 Chapman layer', $
                     ' 4 Z2 - Fit parameter, peak altitude of M2 Chapman layer', $
                     ' 5 H2 - Fit parameter, scale height of M2 Chapman layer', $
                     ' 6 N1 - Fit parameter, peak number density of M1 Chapman layer', $
                     ' 7 Z1 - Fit parameter, peak altitude of M1 Chapman layer', $
                     ' 8 H1 - Fit parameter, scale height density of M1 Chapman layer', $
                     ' 9 Chi Squared - chi-square goodness-of-fit statistic determined by CURVEFIT - see IDL help for definition', $
                     '10 Status - If 0, "good" fit, else "bad" fit, definitions are in fits.pro and fitsfn.pro', $
                     '11 SZA - Solar Zenith Angle, will probably be empty if not MGS profiles', $
                     '12 n0 - number density of lowest altitude data point used in fit', $
                     '13 n1 - number density of highest altitude data point used in fit', $
                     '14 z0 - altitude of lowest altitude data point used in fit', $
                     '15 z1 - altitude of highest altitude data point used in fit', $ 
                     '16 proRMS - RMS of the fit to the profile', $
                     '17 perRMS - RMS of the fit to the profile as a percentage of Nmax' ] 

revflag = KEYWORD_SET(rev)
roundtwoflag = KEYWORD_SET(roundtwo)
IF N_ELEMENTS(tol) EQ 0 THEN tol = 1.0E-3
IF N_ELEMENTS(istart) EQ 0 THEN istart = 0
IF N_ELEMENTS(badfits) EQ 0 THEN badfits = 2
IF N_ELEMENTS(Ain) EQ 0 THEN Ain = [1.0E5,140,10,5E4,110,10]


FOR f = 0,1 DO BEGIN

    bestfits = bestfits - bestfits      ;empty bestfits for second round, just in case
    bestfits[0,*]=INDGEN(nprof)         ;number each profiles
    
;---------------------------------------------------------------------------------------------------------    
;----------------------------------------- FITTING LOOP --------------------------------------------------    
;--------------------------------------------------------------------------------------------------------- 
    
    FOR i=istart,nprof-1 DO BEGIN       ;go through this loop once for each profile

      z1 = altsarr[i,*]        ;altitude scale of electron density profile to be fit
      nonzero = WHERE(z1 NE 0 and z1 lt 200)
      ;************************temporary edit to exclude alts above 200
      z = z1[nonzero]                       ;remove empty elements from the profile
      Nelec = REFORM(eldenarr[i,nonzero])   ;electron density profile to be fit 
      
      if f eq 1 and i eq testi then begin
        window, 0
        plot, nelec, z, /xlog, $
              title = testi, xtitle = 'Electron Density [cm!U-3!N]', ytitle = 'Altitude [km]'
      endif

;-------------------------------------------- Choose Fit Range -------------------------------------------- 

    ;specify which part of the data will be used to run the fitting procedure - will cut off the highest and lowest
    ;part of the electron density profile, specifically the portion above the peak where the number density is less
    ;than ph% and the portion below the peak were the number density is less than pl%

      nmax = MAX(Nelec)
      znmax = z[WHERE(Nelec EQ nmax)]
      znmax = znmax[0]                   ;znmax is a one element array - change it to a number
      bestfits[1,i] = nmax
      zzm = [znmax,znmax]
      nnm = [nmax,nmax]
      
      CASE f of
        0: BEGIN  ;for first round, use the input values
           npl = pl*nmax 
           nph = ph*nmax
           END
        1: BEGIN  ;for the second round, use 50% of the M1 value at the low end
           npl = r2pl*bestfits1[6,i]
           nph = ph*nmax
           END
        ELSE: STOP
      ENDCASE

      if i eq testi then print, npl, nph
    
      zh = z[WHERE(z GE znmax)]          ;split profile into above and below M2 peak
      zl = z[WHERE(z LT znmax)]
      nh = nelec[WHERE(z GE znmax)]
      nl = nelec[WHERE(z LT znmax)]
    
      zfocush = zh[WHERE(nh GE nph)]      ;cut off data at specified percentages to fit to only a section of the data
      zfocusl = zl[WHERE(nl GE npl)]
      focush = nh[WHERE(nh GE nph)]       ;on second round, what is bestfits[6,i] is bad and where gives -1??
      focusl = nl[WHERE(nl GE npl)]
    
      focusz = [zfocusl,zfocush]         ;splice back together the parts that will be fit
      focusn = [focusl,focush]
      
      focus_sort = SORT(focusz)
      focusz = focusz[focus_sort]
      focusn = focusn[focus_sort]

      IF revflag THEN BEGIN
        focusn = reverse(focusn)
        focusz = reverse(focusz)
      ENDIF
      ;if revflag set, switches arrays to descending order to match the original fit procedure
      ;these lines (REVERSE) not in the first fit to the model data  
      ;if revflag off, there is a higher percentage of good fits, so don't recommend setting /REV
      ;actually this is not generally true?
        
      IF N_ELEMENTS(focusn) LE 6 THEN BREAK ;if the relevant portion of the profile is too small, don't fit it
        
      bestfits[14:15,i] = [MIN(focusz,izmin), MAX(focusz,izmax)]
      bestfits[12:13,i]=[focusn[izmin],focusn[izmax]] ;the points at which the fit is cut

      if i eq testi then begin
          help, focusz
          help, focusn
          print, focusz[izmin], focusz[izmax]
          print, focusn[izmin], focusn[izmax]
          print, ''
      endif    

;---------------------------------------------- Do the Fit --------------------------------------------- 
    ;use curvefit to find the best fit parameters
      ;A=[1.0E5,140,10,5E4,110,10]           ;need to give curvefit an initial guess of the fit parameters
;      A=[1.0E5,120,10,5E4,100,10]           ;need to give curvefit an initial guess of the fit parameters
      A = Ain
      ;Ain=A                                 ;curvefit will change the values in A so save a copy as Ain
      Alabel=['N2','Z2','H2','N1','Z1','H1']
    
      WEIGHTS=FLTARR(N_ELEMENTS(focusn))+1
      yfit = CURVEFIT(focusz, focusn, WEIGHTS, A, SIGMA, CHISQ=chi, STATUS=stat, $
          FUNCTION_NAME='chapman', ITMAX=100, YERROR = yerr, /NODERIVATIVE, TOL = tol)
    
      ;sigmas[*,i] = sigma
      ;yerrs[i] = yerr
      bestfits[10,i] = stat       ;save status determined by curvefit
      bestfits[3:8,i]=A
      bestfits[9,i]=chi
      bestfits[2,i]=Chapmanfn(bestfits[4,i],A)               ;fitted Nmax - electron density at M2 altitude

      proRMS = RMS(Chapmanfn(focusz,A), focusn)              ; RMS of profile fit - measure of goodness of fit
      perRMS = proRMS/nmax
      bestfits[16:17,i] = [proRMS, perRMS]    
      
      if i eq testi then begin
          print, i
          setwin
          window,1+f
          plot, focusn, focusz, yrange = [0,250], /xlog, xrange = [1E3,1E6]
          oplot, Chapmanfn(focusz,A), focusz, color = red
          oplot, [focusn[izmin],focusn[izmax]],[focusz[izmin],focusz[izmax]], psym = 2
          print, 'Nmax = '+STRING(nmax)
          print, 'proRMS = '+STRING(prorms)
          print, 'perRMS = '+STRING(perrms)
          print, 'A', bestfits[3:8,i]
          print, ''
          ;stop
      endif
    
    ;...............................Plot the Profile to Check....................................
    IF KEYWORD_SET(showplots) THEN BEGIN
      IF (f eq 1 AND roundtwoflag) OR (f EQ 0 AND NOT roundtwoflag) THEN BEGIN
    
      win = 2
      SET_PLOT, 'win'
      WINDOW, win, xsize = 600, ysize = 400
      LOADCT, 39, /SILENT
      DEVICE, DECOMPOSED=0
      
      PLOT, Nelec, z, xrange = [1E3,1E6], yrange = [0,250], /xlog, $
        xtitle = 'Electron Density [cm!U-3!N]', ytitle = 'Altitude [km]', title = 'Profile '+string(i, format = '(I4)')
        OPLOT, Chapmanfn1(z,bestfits[3:5,i]), z, color=90   ;M2 Chapman layer
        OPLOT, Chapmanfn1(z,bestfits[6:8,i]), z, color=150  ;M1 Chapman layer
        OPLOT, Chapmanfn(z,bestfits[3:8,i]),  z, color=254  ;sum of the two layers
        OPLOT, bestfits[12:13,i], bestfits[14:15,i], psym=2, color=150  ;the two points at which the fit was cut
        ox = 0.9E5
        oy = 240
        oyd = 15
        XYOUTS, ox, oy, 'Model Profile', charsize = 1.5
        XYOUTS, ox, oy-oyd, 'M2 Chapman fit', charsize = 1.5, color = 90
        XYOUTS, ox, oy-2*oyd, 'M1 Chapman fit', charsize = 1.5, color = 150
        XYOUTS, ox, oy-3*oyd, 'Profile fit', charsize = 1.5, color = 254
        formats = ['(E9.2)', '(F9.2)', '(F9.2)', '(E9.2)', '(F9.2)', '(F9.2)']
        FOR j = 0,5 DO XYOUTS, ox, oy-(4+j)*oyd, bestfits_label[j+3]+' = '+string(bestfits[j+3,i], f = formats[j]), $
            charsize = 1.5
        IF N_ELEMENTS(zenarray) NE 0 THEN XYOUTS, ox, oy-11*oyd, 'SZA = '+string(zenarray[i], F='(F5.2)'), charsize = 1.5        
      ;PRINT, 'Profile '+string(i, format = '(I4)')
    
      IF i EQ psind THEN BEGIN
        SET_PLOT, 'PS'
        filename = 'C:\Users\kfallows\IDLWorkspaces\ExampleModelFit'+STRING(psind,F='(I04)')+'.ps'
        PRINT, 'Creating Example Fit Plot: '+filename
        DEVICE, FILENAME = filename, /portrait, /color
        !P.MULTI = 0
        PLOT, Nelec, z, xrange = [0,1.5E5], yrange = [0,250], $
          xtitle = 'Electron Density [cm!U-3!N]', ytitle = 'Altitude [km]', title = 'Profile '+string(i, format = '(I4)'), $
          xthick = 3, ythick = 3, thick = 3, charthick = 3
          OPLOT, Chapmanfn1(z,bestfits[3:5,i]), z, color=90, thick = 3    ;M2 Chapman layer
          OPLOT, Chapmanfn1(z,bestfits[6:8,i]), z, color=150, thick = 3  ;M1 Chapman layer
          OPLOT, Chapmanfn(z,bestfits[3:8,i]),  z, color=254, thick = 3  ;sum of the two layers
          OPLOT, bestfits[12:13,i], bestfits[14:15,i], psym=2, color=150, thick = 3  ;the two points at which the fit was cut
          ox = 0.9E5
          oy = 230
          oyd = 15
          XYOUTS, ox, oy, 'Model Profile', charsize = 1.5, charthick = 3
          XYOUTS, ox, oy-oyd, 'M2 Chapman fit', charsize = 1.5, color = 90, charthick = 3
          XYOUTS, ox, oy-2*oyd, 'M1 Chapman fit', charsize = 1.5, color = 150, charthick = 3
          XYOUTS, ox, oy-3*oyd, 'Profile fit', charsize = 1.5, color = 254, charthick = 3
          FOR j = 0,5 DO XYOUTS, ox, oy-(4+j)*oyd, bestfits_label[j+3]+' = '+string(bestfits[j+3,i], format = '(F9.2)'), charthick = 3
          IF N_ELEMENTS(zenarray) NE 0 THEN XYOUTS, ox, oy-11*oyd, 'SZA = '+string(zenarray[i], F='(F5.2)'), charsize = 1.5
          DEVICE, /CLOSE
          set_plot, 'win'
      ENDIF
    
      ;WAIT, 3
      ;stop
      ENDIF
    ENDIF   ;IF plots EQ 'yes' THEN BEGIN

      if i eq testi then begin
        SET_PLOT, 'PS'
        filename = 'C:\Users\kfallows\IDLWorkspaces\TestiModelFit.ps'
        PRINT, 'Creating Example Fit Plot: '+filename
        DEVICE, FILENAME = filename, /color, /portrait
        !P.MULTI = 0
        PLOT, Nelec, z, xrange = [0,1.5E5], yrange = [0,250], $
          xtitle = 'Electron Density [cm!U-3!N]', ytitle = 'Altitude [km]', title = 'Profile '+string(i, format = '(I4)'), $
          xthick = 3, ythick = 3, thick = 3, charthick = 3
          OPLOT, Chapmanfn1(z,bestfits[3:5,i]), z, color=90, thick = 3    ;M2 Chapman layer
          OPLOT, Chapmanfn1(z,bestfits[6:8,i]), z, color=150, thick = 3  ;M1 Chapman layer
          OPLOT, Chapmanfn(z,bestfits[3:8,i]),  z, color=254, thick = 3  ;sum of the two layers
          OPLOT, bestfits[12:13,i], bestfits[14:15,i], psym=2, color=150, thick = 3  ;the two points at which the fit was cut
          ox = 0.9E5
          oy = 350
          oyd = 15
          XYOUTS, ox, oy, 'Model Profile', charsize = 1.5, charthick = 3
          XYOUTS, ox, oy-oyd, 'M2 Chapman fit', charsize = 1.5, color = 90, charthick = 3
          XYOUTS, ox, oy-2*oyd, 'M1 Chapman fit', charsize = 1.5, color = 150, charthick = 3
          XYOUTS, ox, oy-3*oyd, 'Profile fit', charsize = 1.5, color = 254, charthick = 3
          FOR j = 0,5 DO XYOUTS, ox, oy-(4+j)*oyd, bestfits_label[j+3]+' = '+string(bestfits[j+3,i], format = '(F9.2)'), charthick = 3
          IF N_ELEMENTS(zenarray) NE 0 THEN XYOUTS, ox, oy-11*oyd, 'SZA = '+string(zenarray[i], F='(F5.2)'), charsize = 1.5
          DEVICE, /CLOSE
          set_plot, 'win'
        endif


    ;if f eq 1 then stop
    ENDFOR  ;FOR i=0,nprof-1 DO BEGIN
    
;----------------------------------------------------------------------------------------------------------     
;-------------------------------------------- END FITTING LOOP -------------------------------------------- 
;---------------------------------------------------------------------------------------------------------- 

;---------------------------------------------Bad Fit Definitions------------------------------------------ 
    
    ns = 7
    Status = INTARR(ns,nprof)
    Status[0,*] = bestfits[10,*]
    gi = DBLARR(nprof)

    CASE badfits OF
    
        0: BEGIN 
              
;            PRINT, 'Flagging bad fits'

            ;define parameters for the histogram and Gaussian fit
            typicalvalues = Ain
            maxin = typicalvalues*10
            minin = -maxin
            b=[2.E3,1,0.5,2E3,2,1]    ;binsize
            nbins = (maxin-minin)/b
            
            xlabel=['M2 Peak Electron Density','M2 Peak Altitude','M2 Scale Height', $
                    'M1 Peak Electron Density','M1 Peak Altitude','M1 Scale Height']
            param = ['N2','Z2','H2','N1','Z1','H1']

            ;specified early on in Chap fit analysis
            forcecutlow = [3.0E4, 125, 4, 0, 95, 0]
            forcecuthigh = [1.4E5, 160, 19, 7.0E4, 120, 12]
            
            ;goodrms = WHERE(bestfits[16,*] LE 0.05)
            ;bestfits = bestfits[*,goodrms]

            SETWIN
            !p.charsize = 1.2            
            window,0
            cghistoplot, bestfits[17,*], xrange = [-1,2], mininput = 0, maxinput = 2, binsize = 0.01
            
            window,1

            FOR p=0,nparam-1 DO BEGIN
                
                ;PRINT, ''
                ;PRINT, 'Parameter: '+param[p]
                
                HISTOBIN, bestfits[3+p,*], INDGEN(nbins[p]+1)*b[p] + minin[p], b[p], maxin[p], minin[p], $
                          avgbin, stdevbin, xmid, binelem    
                Gfit = GAUSSFIT(xmid, binelem, GA, CHISQ = chisq, NTERMS = 3, SIGMA = sig)
                
                PLOT, xmid, binelem, xtitle = xlabel[p], ytitle = 'Frequency', xrange = [-10,10]*GA[2]+GA[1]
                    OPLOT, xmid, Gfit, color = red
                
                    nsig = 5.0
                    VLINE, GA[1], color = red
                    VLINE, GA[1]-nsig*GA[2], color = blue
                    VLINE, GA[1]+nsig*GA[2], color = blue
                    VLINE, forcecutlow[p], color = green
                    VLINE, forcecuthigh[p], color = green
                    
                    AL_LEGEND, ['Histogram', 'Gauss Fit', 'Force Cut', +STRING(nsig,F='(I1)')+' Sigma Cut'], $
                        linestyle=0, colors = [black, red, green, blue], box = 0
                
                
                mostpop = MAX(binelem, imostpop)    ;most populated bin of parameter
                
                PRINT, '  max fit value  = '+STRING(maxin[p])
                PRINT, '  min fit value  = '+STRING(minin[p])
                PRINT, '  most populated = '+STRING(xmid[imostpop])
                PRINT, '  Gauss center   = '+STRING(GA[1])
                PRINT, '  # in most pop. = '+STRING(mostpop)
                PRINT, '  Gauss height   = '+STRING(GA[0])
                PRINT, '  Gauss st. dev. = '+STRING(GA[2])
                PRINT, ''
                
                FOR i=0,nprof-1 DO BEGIN
                
                  IF p EQ 5 THEN $
                  IF bestfits[3+p,i] LT GA[1]-nsig*GA[2] OR bestfits[3+p,i] GT GA[1]+nsig*GA[2] THEN Status[1+p,i]=3+p
                  IF bestfits[3+p,i] LT 0.0 OR bestfits[3+p,i] GT GA[1]+nsig*GA[2] THEN Status[1+p,i]=3+p
                  ;if a fit is identified as bad by above block of if statements, then add the flag "3"
                  ; to the status column of bestfits, but don't replace a status of 1 or 2 from curvefit
                  IF MAX(status[*,i]) GE 3 THEN $
                  IF bestfits[10,i] EQ 0 THEN bestfits[10,i] = 3

                ;If the N1 was bad the first round, don't trust the fit the second time around and replace with the
                ;first fit, which will be flagged as a bad fit.  
                
                IF f EQ 1 THEN $
                IF status1[4,i] NE 0 THEN BEGIN
                    bestfits[*,i] = bestfits1[*,i]
                    status[*,i] = status1[*,i]
                ENDIF
                ENDFOR
                
                IF f EQ 0 THEN status1 = status
                
                wait,2
                ;stop
            ENDFOR
    
            END   ;CASE badfits OF 0
            
        1: BEGIN    ;CASE badfits
        
          PRINT, 'BADFITS = 1'
        
            FOR i=0,nprof-1 DO BEGIN
              ;Chopping off each parameter within a certian range, based largely on the 
              ;distribution of values in a histogram of each parameter
              ;NOTE these are very specific to MGS set. Not valid at lower SZA, for example
              IF bestfits[3,i] LT 3.0E4 OR bestfits[3,i] GT 1.4E5 THEN Status[1,i]=3    
              IF bestfits[4,i] LT   125 OR bestfits[4,i] GT   160 THEN Status[2,i]=4 
              IF bestfits[5,i] LT     4 OR bestfits[5,i] GT    19 THEN Status[3,i]=5 
              IF bestfits[6,i] LT     0 OR bestfits[6,i] GT 7.0E4 THEN Status[4,i]=6
              IF bestfits[7,i] LT    95 OR bestfits[7,i] GT   120 THEN Status[5,i]=7
              IF bestfits[8,i] LT     0 OR bestfits[8,i] GT    12 THEN Status[6,i]=8

              ;if a fit is identified as bad by above block of if statements, then add the flag "3"
              ; to the status column of bestfits, but don't replace a status of 1 or 2 from curvefit
              IF MAX(status[*,i]) GE 3 THEN BEGIN 
                  IF bestfits[10,i] EQ 0 THEN bestfits[10,i] = 3
              ENDIF ;ELSE PRINT, MAX(bestfits[10,i]), bestfits[3:8,i]
 
              if i eq testi then print, status[*,i], bestfits[6,i], bestfits[10,i]
               
            ENDFOR

            IF f EQ 0 THEN status1 = status
              
            END   ;CASE badfit of 1
            
        2: BEGIN    ;CASE badfits

            FOR i=0,nprof-1 DO BEGIN

              ;Alternate set of rules for eliminating bad fits
              ;   - any parameters are negative
              ;   - the M2 altitude is above the top of focusz
              ;   - the M2 scale height is larger than 20 km
              ;   - the M1 altitude is above the M2 altitude, or below the bottom of focusz
              ;   - the M1 scale height is larger than 12 km
              
              checkforneg = WHERE(bestfits[3:8,i] LE 0, negcount)
              IF negcount      GT 0              THEN status[1,i] = 3  
              IF bestfits[4,i] GT bestfits[15,i] THEN status[2,i] = 4                          
              IF bestfits[5,i] GT 20             THEN status[3,i] = 5
              IF bestfits[6,i] GT bestfits[3,i]  THEN status[4,i] = 6
              IF bestfits[7,i] LT bestfits[14,i] OR bestfits[7,i] GT bestfits[4,i] THEN status[5,i] = 7
              IF bestfits[8,i] GT 12             THEN status[6,i] = 8

              ;if a fit is identified as bad by above block of if statements, then add the flag "3"
              ; to the status column of bestfits, but don't replace a status of 1 or 2 from curvefit
              IF MAX(status[*,i]) GE 3 THEN $
              IF bestfits[10,i] EQ 0 THEN bestfits[10,i] = 3
 
              ;THIS TAKES OUT 200 PROFILES. NOT SURE THAT I ACTUALLY WANT TO DO THIS.
              ;if in the first round, the M1 layer was located incorrectly, keep it as bad
              ;the lower boundary of the fit region cannot be identified reasonably if the 
              ;M1 fit parameters are too crazy this first time around
              ;(trying to avoid WHERE(nl GE npl) = -1)
              IF f EQ 1 THEN BEGIN
                  IF status1[5,i] NE 0 THEN BEGIN
                      bestfits[*,i] = bestfits1[*,i]
                      status[*,i] = status1[*,i]
                      ;need to replace status explicitly, because the CURVEFIT status 1 and 2
                      ;will not get automatically set by going through the badfits loop again
                  ENDIF
              ENDIF
              
              IF f EQ 1 THEN BEGIN 
              if i eq testi then begin

              print, bestfits[3:8,testi]
              print, bestfits1[3:8,testi]
              endif
              endif

              if i eq testi then print, status[*,i], bestfits[6,i], bestfits[10,i]
              
            ENDFOR  ;i=0,nprof-1

            IF f EQ 0 THEN status1 = status
              
            END

          
        ELSE:  PRINT, 'No bad fits flagged'
          
    ENDCASE   ;badfits
        
    FOR i=0,nprof-1 DO BEGIN
      gc=0
      FOR g=0,ns-1 DO BEGIN                 
        IF status[g,i] NE 0 THEN gc=gc+1
        gi[i]=gc
      ENDFOR
    ENDFOR
    
    good=WHERE(gi EQ 0)
    bad=WHERE(gi NE 0)
    ng=N_ELEMENTS(good)
    nb=N_ELEMENTS(bad)
    
    nbi=intarr(8)-1                                             ;bi is a 2D array containing the indices in bestfits of 
    bi=intarr(8,nb)-1                                           ;all of the bad fits - one column for each bad fit definition
    bi0=WHERE(status[0,*] EQ 1)
    nbi[0]=N_ELEMENTS(bi0)                                      ;nbi is an array containing the number of bad fits of each type 
    FOR i=1,7 DO BEGIN                                          ;- one element for each bad fit definition
    bi[0,0:nbi[0]-1]=bi0
      bii=WHERE(status[i-1,*] EQ i+1)
      nbii=N_ELEMENTS(bii)
      bi[i,0:nbii-1]=bii                                                                    
      nbi[i]=nbii                                                                           
    ENDFOR
    ;ENDIF

              ;i = 3481
              ;print, status[*,i], bestfits[6,i], bestfits[10,i]

    IF NOT roundtwoflag THEN BREAK    ;only do the second round if specified
    bestfits1 = bestfits              ;save copy from first round
    good1 = good
    ;IF f EQ 0 THEN PRINT, 'Beginning Second Round'
    
ENDFOR    ;f
              ;i = 3481
              ;print, status[*,i], bestfits[6,i], bestfits[10,i]

IF N_ELEMENTS(savefile) NE 0 THEN SAVE, bestfits, good, bestfits1, good1, bestfits_label, bestfits_comments, $
                                        badfits, ph, pl, r2pl, filename = savefile
;stop
;SAVE, bestfits, bestfits_label, bestfits_comments, good, bad, fitsg, fitsb, nbi, bi, jdates, sigmas, yerrs, $
;  filename=path+'Chapman Fitting\Fit Parameters 5600\mgschapfits_gb.sav'


END