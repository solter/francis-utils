##-*- texinfo -*-
##@deftypefn{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} agErlyDef(@var{A},@var{topshifts})
##@deftypefnx{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} agErlyDef(@var{A},@var{topshifts},@var{toplt})
##
##This introduces 2 bulges (1 from top and 1 from bottom) which sweep until they meet.@*
##Then a Schur decomposition is performed to create spikes replacing the bulges.@*
##
##Inputs:@*
##  @var{A} - A hessenberg matrix. Results will be meaningless if not hessenberg.@*
##  @var{topshifts}, - Lists of shifts to introduce in the bulges respectively.
##    Should be less than ~6 to avoid shift blurring.@*
##  @var{toplt} - Optional argument. If it is true, this function
##    will plot each step of the bulge creation and chasing.@*
##  @var{toprt} - Optional argument. If it is true, this function
##    will create a directory impStepPlts save the plots to it.@*
##
##Outputs:@*
##  @var{H} - The matrix A after running 2 bulges into it's center and
##    creating the spikes@*
## @end deftypefn
function [H, pltNum] = agErlyDef(A, Shifts, toplt = false, toprt = false)
  H=A;

  sl = length(Shifts);
  pl = .2;
  
  if(toplt)
    mkdir('impStepPlts');
    pltNum = 0;
  else
    toprt = false;
  end#if

  #put shifts into H
  #add bulge to the top
  for i=1:sl
    endIdx = 1+i;
    
    #create house vector
    v = H(1:endIdx, 1);
    v(1) -= Shifts(i);
    v(1) += sgn(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);#normalize house vector
    #apply householder transformation to the right bits
    H(1:endIdx,:) -= v*((2*v')*H(1:endIdx,:));
    H(:,1:endIdx) -= (H(:,1:endIdx)*(2*v))*v';#many zero mulitplies

    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(pl);
      end#if
    end#if
  endfor


  #move bulge 1-column over per i
  stIdx = 1;
  endIdx = 1 + sl;
  do
    ++stIdx;
    ++endIdx;
    
    #move top bulge
    #create house vector
    v = H(stIdx:endIdx, stIdx-1);
    v(1) += sgn(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);#normalize house vector
    #apply householder transformation to the right bits
    H(stIdx:endIdx,:) -= v*((2*v')*H(stIdx:endIdx,:));
    H(1:endIdx+1,stIdx:endIdx) -= (H(1:endIdx+1,stIdx:endIdx)*(2*v))*v';
    H(stIdx+1:endIdx,stIdx-1) = 0;#zeros everything out for exactness
   
    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(pl);
      endif
    end#if

  until(endIdx + 1 >= length(H))#if it touches the bottom 


  #create spikes
  spSt = stIdx+1;#index of block
  [spRot,~] = schur(H(spSt:end,spSt:end));
  H(:,spSt:end) = H(:,spSt:end)*spRot;
  H(spSt:end,(spSt-1):end) = spRot'*H(spSt:end,(spSt-1):end);

  if(toplt)
    pltMat(H);
    if(toprt)
      print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
    else
      pause(pl);
    endif
  end#if

  if(toprt)
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
