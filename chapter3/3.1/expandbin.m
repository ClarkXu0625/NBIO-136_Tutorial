function [new_vector] = expandbin(old_vector, old_dt, new_dt)
% The function will take as inputs an initial vector (old_vector), an 
% initial bin-width (old_dt), and a final bin-width (new_dt). 
% It will return an output vector of smaller size than the initial vector, 
% with the size reduced by a factor equal to the ratio of the bin-widths.

    ratio = round(new_dt/old_dt);   % times new vector greater than old

    Nold = length(old_vector);  % length of old vector
    
    Nnew = round(Nold/ratio);  % length of new vector

    new_vector = zeros(1, Nnew);   

    for i=1:Nnew
        new_vector(i) = sum(old_vector((i-1)*ratio+1 : i*ratio))/ratio;        
    end
end