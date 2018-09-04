function Iout = f_StackProjections(SEP,Type)

    if  strcmp(Type,'Maximum')
        sizeX = size((SEP(1).data), 1);
        sizeY = size((SEP(1).data), 2);
        Iout = SEP(1).data;
        h = waitbar(0,strcat ('Please wait computing maximumprojection'));
        for a=1:(size(SEP,2))-1
            frameAdd = SEP(a+1).data;
            for i =1 : sizeX
                for j=1: sizeY
                    if frameAdd(i, j) > Iout(i, j)
                        Iout (i,j) = frameAdd(i,j);
                    end
                end
            end
            waitbar(a/(size(SEP,2)))
        end
        Iout  = uint16(Iout);
        close(h)

    elseif strcmp(Type,'Sum')
        Iout = (zeros(size(SEP(1).data,1),size(SEP(1).data,2)));
        for i=1:size(SEP,2)
            Iout = imadd(Iout,double(SEP(i).data));
        end
        
    elseif strcmp(Type,'Mean')
        Iout = zeros(size(SEP(1).data,1),size(SEP(1).data,2));
        for i=1:size(SEP,2)
            Iout = imadd(Iout,double(SEP(i).data));
        end
        Iout = uint16(Iout/ size(SEP,2));
    end
    
end

