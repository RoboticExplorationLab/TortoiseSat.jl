function DatetoJD(Y,M,D,H,M,S)
#this function converts traditional units of measurement into Julian Day UTC

JDN=(1461 × (Y + 4800 + (M − 14)/12))/4 +(367 × (M − 2 − 12 × ((M − 14)/12)))/12 − (3 × ((Y + 4900 + (M - 14)/12)/100))/4 + D − 32075

end
