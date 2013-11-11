%H must be hessenberg, shifts should be a vector
%OUTPUTS: H - 1-sweep after H
function [H] = impStep(A, shifts, toplt)
  H=A;
  sl = length(shifts);
  bs = sl+1;%bulge size
  
  if(toplt)
    mkdir('impStepPlts');
    pltNum = 0;
  end%if

  %put shifts into H
  for i=1:sl
    %fast version
    endIdx = 1+i;
    %create house vector
    v = H(1:endIdx, 1);
    v(1) += sign(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);%normalize house vector
    %apply householder transformation to the right bits
    H(1:endIdx,:) -= v*((2*v')*H(1:endIdx,:));
    H(:,1:endIdx) -= (H(:,1:endIdx)*(2*v))*v';
    if(toplt)
      pltMat(H);
      print(sprintf('impStepPlts/impstep%03d.png',++pltNum));
    end%if
  end%for

  %move bulge 1-column over per i
  for i=1:length(H)-2
    endIdx = min(i + bs,length(H));
    %create house vector
    v = H(i+1:endIdx, i);
    v(1) += sign(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);%normalize house vector
    
    Q = eye(length(v),length(v)) - 2*v*v';
    if(norm(Q*Q - eye(length(v),length(v)),'inf') > 1e-8)
      printf('dafuq at i = %02d\n', i);
    end

    %apply householder transformation to the right bits
    H(i+1:endIdx,:) -= v*((2*v')*H(i+1:endIdx,:));
    H(1:min(endIdx+1,length(H)),i+1:endIdx) -= (H(1:min(endIdx+1,length(H)),i+1:endIdx)*(2*v))*v';
    %zero out appropriate bits
    H(i+2:min(i+bs,length(H)),i) = 0;%clean up to ensure 0
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
  end%for

  if(toplt)
    close all;
  end%if
end%function

function [] = pltMat(H)
    imagesc(log(abs(H)) + .01);
end%function
