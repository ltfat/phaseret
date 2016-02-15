function E = magnitudeerr(target,reconstructed)

E = norm(abs(target)-abs(reconstructed),'fro')/norm(abs(target),'fro');