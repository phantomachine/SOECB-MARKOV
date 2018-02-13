function [x]=makezero(x)

% Code written by Richard Dennis
% return near zero elements (<1e-9) in the matrix x as zeros.

eta = 1e-9;
i = 1;
j = 1; 
nr = rows(x);
nc = cols(x);

while i <= nr;
  j = 1;
  while j <= nc;
    if abs(x(i,j)) < eta;
      x(i,j) = 0;
    end;
    j = j + 1;
  end;
  i = i + 1;
end;
