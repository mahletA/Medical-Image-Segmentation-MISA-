function [e] = SD(Ifixed, img)
     e=sum((img(:)-Ifixed(:)).^2); %/numel(img);
end