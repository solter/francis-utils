function [val] = cbWrap(shifts, A)

  n = length(shifts)/2;

  val = zeros(1,n);
  H = colAndBust(A, shifts(1:n), shifts(n:end));

  colIdx = floor(length(A)/2) - (n + 1);
  rowIdx = colIdx + 2*n + 3;
  val(1:n) = H(rowIdx - n-1:rowIdx-1,colIdx)';
  val(2:n) = H(rowIdx,colIdx + 1: colIdx + n + 1);

endfunction
