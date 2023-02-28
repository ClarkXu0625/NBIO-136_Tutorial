function [new_vector] = expandbin(old_vector, old_dt, new_dt)

    ratio = round(new_dt/old_dt);   % times new vector greater than old

    Nold = length(old_vector);  % length of old vector
    
    Nnew = round(Nold/ratio);  % length of new vector

    new_vector = zeros(1, Nnew);   

    for i=1:Nnew
        new_vector(i) = sum(old_vector((i-1)*ratio+1 : i*ratio))/ratio;        
    end


end