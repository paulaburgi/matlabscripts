function [noise,covd]=make_corr_noise(x,y,Var,noisescale,numdatas,weights,pixelsize)

if(1)
  np=length(x(:));
  covd=diag(ones(np,1)*Var);
  if(noisescale>0)
      for i=1:np
          
          tmp=sqrt((x(i)-x(:)).^2+(y(:)-y(i)).^2);
          covd(i,:)=Var*10.^(-tmp/noisescale);
          
      end
  end
  
else
  np=length(x);
  
  h=waitbar(0,'fast variance');
  for i=1:np
    bxx      = -weights(i)/2:weights(i)/2;
    bxy      = -weights(i)/2:weights(i)/2;
    
    nx       = length(bxx);
    tx       = 0:nx-1;
    
    clear sums dists mults
    for j=1:nx
      for k=1:nx
	dists(k,:) = sqrt((tx-tx(k)).^2+tx(j)^2);
      end
      covs     = Var*10.^(-dists*pixelsize/noisescale);
      sums(j)  = sum(covs(:));
      mults(j) = [nx-j+1]*2;
    end
    mults(1)=nx;
    covd(i,i)=sums*mults'/nx^4;
    
    for l=i+1:np
      dists = sqrt((x(i)-x(l))^2+(y(i)-y(l))^2);
      covd(i,l)=Var*10.^(-dists/noisescale);
      covd(l,i)=covd(i,l);
    end
    waitbar(i/np,h);
  end
  close(h);
end

noise=corr_noise(covd,numdatas);
