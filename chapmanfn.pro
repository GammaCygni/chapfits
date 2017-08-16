;Katy Fallows
;Boston University

function Chapmanfn, z, A 

;----------------------------------------------------------------
;This is the same equation as chapman.pro, used to fit the 
;MGS electron density profiles.
;
;INPUT
; z - an altitude scale
; A - 6 element array containing Nm2, zm2, H2, Nm1, zm1, H1
;   A[0] is the electron density at the peak of the Chapman layer
;   A[1] is the altitude of the peak of the Chapman layer
;   A[2] is the scale height of the Chapman layer
;   A[3]-A[5] are the same values for a second layer 
;OUTPUT
; N - electron density profile
;----------------------------------------------------------------

N2=A[0]*exp(.5*((1-(z-A[1])/A[2])-exp(-(z-A[1])/A[2])))
N1=A[3]*exp(.5*((1-(z-A[4])/A[5])-exp(-(z-A[4])/A[5])))
N=N2+N1

return, N

end