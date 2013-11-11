%This introduces 2 bulges (1 from top and 1 from bottom)
%which sweep until they meet.
%Then a Schur decomposition is performed to bust it open
%and create spikes. Ideally this will reorder it to sort
% the spikes lowest to hightest
%
%H must be hessenberg, shifts should be a vector
%OUTPUTS: The busted open shifts
function [H] = colAndBust(A, shifts, toplt)
  H=A;
  sl = length(shifts);
  
  if(toplt)
    mkdir('impStepPlts');
    pltNum = 0;
  end%if

  %put shifts into H
  for i=1:sl
    endIdx = 1+i;
    
    %add bulge to the top
    %create house vector
    v = H(1:endIdx, 1);
    v(1) -= shifts(i);
    v(1) += sign(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);%normalize house vector
    %apply householder transformation to the right bits
    H(1:endIdx,:) -= v*((2*v')*H(1:endIdx,:));
    H(:,1:endIdx) -= (H(:,1:endIdx)*(2*v))*v';%many zero mulitplies

    %add bulge to the bottom
    v = H(end,(end-endIdx):end);
    v(end) -= shifts(i);
    v(end) += sign(v(end))*sqrt(v*v');
    v /= sqrt(v*v');%normalize house vector
    H((end-endIdx):end,:) -= v'*((2*v)*H((end-endIdx):end,:));%many zero multiplies
    H(:,(end-endIdx):end) -= (H(:,(end-endIdx):end)*(2*v'))*v;
    
    if(toplt)
      pltMat(H);
      print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
    end%if
  end%for

  %move bulge 1-column over per i
  n = length(H);
  i=1;
  do
    stIdx = i+1;
    endIdx = i+1 + sl;
    
    %move top bulge
    %create house vector
    v = H(stIdx:endIdx, i);
    v(1) += sign(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);%normalize house vector
    %apply householder transformation to the right bits
    H(stIdx:endIdx,:) -= v*((2*v')*H(stIdx:endIdx,:));
    H(1:endIdx+1,stIdx:endIdx) -= (H(1:endIdx+1,stIdx:endIdx)*(2*v))*v';
    H(stIdx+1:endIdx,i) = 0;%zeros everything out for exactness
    
    %move bottom bulge
    endIdx = n-i;
    stIdx = n-i-sl;
    v = H(endIdx + 1,stIdx:endIdx);
    v(end) += sign(v(end))*sqrt(v*v');
    v /= sqrt(v*v');%normalize house vector
    H(stIdx:endIdx,stIdx - 1:n) -= v'*((2*v)*H(stIdx:endIdx,stIdx-1:n));
    H(:,stIdx:endIdx) -= (H(:,stIdx:endIdx)*(2*v'))*v;
    H(endIdx + 1,stIdx:endIdx-1) = 0;%zeros stuff out for exactness
    
    i++;
    if(toplt)
      pltMat(H);
      print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
    end%if
    %create spikes
%    if(0)%(i + 1 + length(shifts) < length(H))
%    blkIdx = (i+1):(i+1+length(shifts));
%    subH = H(blkIdx,blkIdx);
%    [r1,~] = schur(subH);
%    bigr = eye(size(H));
%    bigr(blkIdx,blkIdx) = r1;
%    pltMat(bigr'*H*bigr);
%    SpikeFrame(i) = getframe;
%    end%if
  until(i+1+sl >= n-i-sl)

  if(toplt)
    close all;
  end%if
end%function

function [] = pltMat(H)
    imagesc(log(abs(H)) + .01);
end%function
