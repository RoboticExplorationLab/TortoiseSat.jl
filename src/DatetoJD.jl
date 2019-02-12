function DatetoJD(Y,M,D,Hour,Min,Sec)
#this function converts traditional units of measurement into Julian Day UTC

JD=367*Y-(7*(Y+5001+(M-9)/7))/4+(275*M)/9+D+1729777
JD += (Hour-12)/24+(Min)/1440+(Sec)/86400
end
