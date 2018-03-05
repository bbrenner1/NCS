function [ out ] = fftshift2( x )
%fftshift2 is a recreation of the matlab function fftshift
%   fftshift2 circularly shifts each row down by half as many rows as are in the matrix, then
%   shifts each column right by half as many columns as are in the matrix.
%   If the row being shifted down is at the bottom, it gets moved up to the
%   first row. Similarly, if the column being shifted right is already the last column, 
%   it will be moved to the first column.
shift =floor(size(x)/2);
down_shift=shift(1);
right_shift=shift(2);

N=size(x);

temp1=zeros(N);
out=zeros(N);
for n=1:1:N(1)
    if n+down_shift <=N(1)
        temp1(n+down_shift,:)=x(n,:);

    else
        temp1((n+down_shift-N(1)),:)=x(n,:);
    end
end

for n=1:1:N(1)
    if n+right_shift <=N(2)
    
        out(:,n+right_shift)=temp1(:,n);
    else
        out(:,n+right_shift-N(2))=temp1(:,n);
    end
end

end
