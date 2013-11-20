##usage: [H] = colAndBust(A,topshifts, botShifts,toplt = false)
##
##This introduces 2 bulges (1 from top and 1 from bottom)
##which sweep until they meet.
##Then a Schur decomposition is performed to bust it open
##and create spikes. Ideally this will reorder it to sort
## the spikes lowest to hightest
##
##H must be hessenberg, shifts should be a vector
##OUTPUTS: The busted open shifts
function [H] = colAndBust(A, topShifts, botShifts, toplt=false)
  H=A;

  tsl = length(topShifts);
  bsl = length(botShifts);
  
  if(toplt)
    mkdir('impStepPlts');
    pltNum = 0;
  end#if

  #put shifts into H
  
  #add bulge to the top
  for i=1:tsl
    endIdx = 1+i;
    
    #create house vector
    v = H(1:endIdx, 1);
    v(1) -= topShifts(i);
    v(1) += sgn(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);#normalize house vector
    #apply householder transformation to the right bits
    H(1:endIdx,:) -= v*((2*v')*H(1:endIdx,:));
    H(:,1:endIdx) -= (H(:,1:endIdx)*(2*v))*v';#many zero mulitplies

    if(toplt)
      pltMat(H);
      print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
    end#if
  endfor

  #add bulge to the bottom
  for i=1:bsl
    endIdx = 1+i;
    v = H(end,(end-endIdx):end);
    v(end) -= botShifts(i);
    v(end) += sgn(v(end))*sqrt(v*v');
    v /= sqrt(v*v');#normalize house vector
    H((end-endIdx):end,:) -= v'*((2*v)*H((end-endIdx):end,:));#many zero multiplies
    H(:,(end-endIdx):end) -= (H(:,(end-endIdx):end)*(2*v'))*v;
    
    if(toplt)
      pltMat(H);
      print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
    end#if
  end#for

  #move bulge 1-column over per i
  n = length(H);
  tstIdx = 1;
  tendIdx = 1 + tsl;
  bendIdx = n;
  bstIdx = n - bsl;
  do
    ++tstIdx;
    ++tendIdx;
    
    #move top bulge
    #create house vector
    v = H(tstIdx:tendIdx, tstIdx-1);
    v(1) += sgn(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);#normalize house vector
    #apply householder transformation to the right bits
    H(tstIdx:tendIdx,:) -= v*((2*v')*H(tstIdx:tendIdx,:));
    H(1:tendIdx+1,tstIdx:tendIdx) -= (H(1:tendIdx+1,tstIdx:tendIdx)*(2*v))*v';
    H(tstIdx+1:tendIdx,tstIdx-1) = 0;#zeros everything out for exactness
   
    --bendIdx;
    --bstIdx;
    if(tendIdx + 1 < bstIdx)#if they won't intersect
      #move bottom bulge
      v = H(bendIdx + 1,bstIdx:bendIdx);
      v(end) += sgn(v(end))*sqrt(v*v');
      v /= sqrt(v*v');#normalize house vector
      H(bstIdx:bendIdx,bstIdx - 1:n) -= v'*((2*v)*H(bstIdx:bendIdx,bstIdx-1:n));
      H(:,bstIdx:bendIdx) -= (H(:,bstIdx:bendIdx)*(2*v'))*v;
      H(bendIdx + 1,bstIdx:bendIdx-1) = 0;#zeros stuff out for exactness
    else
      bendIdx++;
      bstIdx++;;
    end#if

    if(toplt)
      pltMat(H);
      print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
    end#if

  until(tendIdx + 2 >= bstIdx)#if they are touching


  #create spikes
  spSt = tstIdx+1;#index of left block
  spEnd = bendIdx-1;
  [spRot,~] = schur(H(spSt:spEnd,spSt:spEnd));
  H(1:spEnd+1,spSt:spEnd) = H(1:spEnd+1,spSt:spEnd)*spRot;
  H(spSt:spEnd,(spSt-1):n) = spRot'*H(spSt:spEnd,(spSt-1):n);

  if(toplt)
    pltMat(H);
    print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
  end#if


  if(toplt)
    close all;
  end#if

end#function

function [sn] = sgn(v)
  if(v >= 0)
    sn= 1;
  else
    sn= -1;
  endif
endfunction

function [] = pltMat(H)
    imagesc(log(abs(H)) + .01);
    daspect([1 1 1]);
end#function
