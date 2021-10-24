function CompareSystems(Matrix1, Matrix2, Vector1, Vector2)
close all

tol = 1e-3;

figure(1)
Matdif = Matrix2-Matrix1;
[valmat, pos] = max(max(abs(Matdif)))

check = abs(Matdif)> tol;
Matdif = Matdif .* check;
spy(Matdif);

%figure(2)
Vecdif = Vector2-Vector1;

[val, pos] = max(abs(Vecdif))

check = abs(Vecdif)> tol;
Vecdif = Vecdif .* check;



    
%     if abs(Diff(i)<0.5 || isnan(Diff(i)) || abs(VecSerial(i)) < 1e-2 )
%         Diff(i) = 0;
%     end
% end
% [Vecerror2,posi] = max(abs(Diff))