function mash=mask(a)

mask = imread(a);  %Ïû³ı±ßÔµ»Ô¹âÓ°Ïì
% mash(mask==0) = 0;
mash = mask;
for i=11:1:584-11
    for j=11:1:565-11
        if(mask(i-10,j-10)==0||mask(i+10,j+10)==0)
            mash(i,j) = 0;
        end
    end
end