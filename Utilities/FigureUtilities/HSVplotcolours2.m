function M = HSVplotcolours2(dim1, dim2, dim3)
% define plot colours in HSV colour-space over 3 dimensions
% dim1, dim2, and dim3 are the respective sizes of each dimension
% 
% Hue (H) is used for dimension 1, with equally distributed hues around the
% hue circle.
%
% Value (V) is used for dimension 2.
%
% Saturation (S) is used for dimension 3.
%
% The output is a 4-dimensional matrix, with the RGB values in the fourth
% dimension





% Hue
if dim1 == 1
    H = 0;
else
    H = 0:1/(dim1):(1-1/dim1);
end
H = repmat(H',[1,dim2,dim3]);



% Value
maxV = 0.9;
minV = 0.4;
if dim2 == 1
    V = maxV;
else
    V = minV:(maxV-minV)/(dim2-1):maxV;
end
V = repmat(V,[dim1,1,dim3]);



% Saturation
maxS = 1;
minS = 0.3;
if dim3 ==1;
    S = maxS;
else
    S = minS:(maxS-minS)/(dim3-1):maxS;
end
S = permute(S,[1,3,2]);
S = repmat(S,[dim1,dim2,1]);




% convert HSV to RGB
M = zeros(dim1,dim2,3,dim3);
for z = 1:dim3
    M(:,:,:,z) = hsv2rgb(cat(3,H(:,:,z),S(:,:,z),V(:,:,z)));
end
M = permute(M,[1,2,4,3]);

end % eof