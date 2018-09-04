function [Roi,border_coord] = cut_square_on_Image(Ided, I, squareD)
            
           %cut_square_on_Image
%             squareD=8;
            a = (-squareD:squareD)';
            b= zeros(size(a,1),1)+squareD; 
            c = [a b; a -b; b a; -b a];
            border_vector = unique(c, 'rows');
            %cut square
            current_coord = [round(Ided(1,1)),round(Ided(1,2))];
            %Determine the coordinates of the pixels around this pixel
            border_coord = [border_vector(:,1) + current_coord(1),border_vector(:,2) + current_coord(2)];
            border_coord = border_coord(find(border_coord(:,1) >0  & border_coord(:,2) >0 & border_coord(:,2) <= size(I,1) & border_coord(:,1) <= size(I,2)),:);
            %regions inside r
            Roi = I(min(border_coord(:,2)):max(border_coord(:,2)), min(border_coord(:,1)):max(border_coord(:,1)));
end