
function FILENAME = getname(filename,n)

    temp = filename;
    for I=1:n-1 
        temp = [temp,'_'];
    end
    FILENAME = [temp,'.h5'];
end