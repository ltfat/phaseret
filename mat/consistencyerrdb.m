function E = consistencyerrdb(c,projc)

E = 20*log10(norm(c-projc,'fro')/norm(c,'fro'));