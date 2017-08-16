;the same as chapmanfn.pro, but for a single Chapman layer

function Chapmanfn1, z, A           ;A is 3 element array containing Nm, zm, H

N=A[0]*exp(.5*((1-(z-A[1])/A[2])-exp(-(z-A[1])/A[2])))

return, N

end