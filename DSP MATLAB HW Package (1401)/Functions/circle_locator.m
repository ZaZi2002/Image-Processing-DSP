function center_radius = circle_locator(edge_points , image_size , degree_tolerance)
%
% center_radius = circle_locator(edge_points , image_size , degree_tolerance)
%
% This function implements the gradient-pair method for detecting circles
% in an image. The basis of the method is the following fact:
%
%       the diagonally opposite points on the perimeter of a circle
%       have the same gradient orientation but opposite direction,
%       while the orientation coincides with the line connecting them.
%
% Therefore, is this method, we check the gradient direction of all edge 
% point pairs of the image to see whether they satisfy the above
% constraint. In the positive case, their middle point is a candidate for a
% circle center while half the connecting length is a potential radius. The
% output of this function is a 3D matrix contating all possible
% center/radius candidates. In more details, the output is an r*c*250
% matrix where r,c are the size of the input image. The (i,j,k) element of
% the output indicates the radius of a circle for which the point (i,j) of 
% the image was selected as the center for the k'th time. If (i,j) is
% selected less than k times as a center, this element is set to 0 (a circle 
% with 0 radius). In fact, 250 is an upper-bound on the number of times a
% point could be selected as a center.
%
% "edge_points":
% is the 4*n matrix where 'n' is the number of detected edge pixels/points. 
% Each column of this matrix respectively contains the row, column, gradient 
% magnitude, and gradient direction (in degrees) of a pixel passing the
% threshold test.
%
% "image_size":
% is a 1*2 vector of positive integers representing the size of the input
% image.
%
% "degree_tolerance":
% is a positive real number that indicates the tolerance in degrees for
% identifying the points with similar gradient directions. In simple words,
% instead of the points having exactly the same gradient directions, we
% allow for small deviations (e.g., 1 degree).
%
% 
% "center_radius":
% is a 3D matrix contating all possible center/radius candidates. In more 
% details, "center_radius" is an r*c*250 matrix where "image_size = [r,c]".
% The (i,j,k) element of "center_radius" indicates the radius of a circle 
% for which the point (i,j) of the image serves as the center for the k'th 
% time. If (i,j) is selected less than k times as a center, this element is 
% set to 0 (a circle with 0 radius). In fact, 250 is an upper-bound on the 
% number of times a point could be selected as a center.




% defining the output
center_radius               = zeros( image_size(1) , image_size(2) , 250);

% the number of times each point is selected as a center
selection_counter           = sum(sign(center_radius) , 3);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                            You should modify here                              %%%%%%%
%%%%%%%                                                                                %%%%%%%
%%%%%%%                                      _||_                                      %%%%%%%
%%%%%%%                                      \  /                                      %%%%%%%
%%%%%%%                                       \/                                       %%%%%%%
%%%%%%%                                                                                %%%%%%%
                                                                                       %%%%%%%
                                                                                       %%%%%%%
                                                                                       %%%%%%%
% for checking the edge point pairs, we sweep over all points                          %%%%%%%
h   = waitbar(0 , 'Gradient-pair method');
for edge_points_ind = 1 : size( edge_points , 2) - 1                                     
    waitbar( edge_points_ind / (size(edge_points , 2)) )
    direction1 = edge_points(4,edge_points_ind); % direction of the first pixel
    for edge_points_ind2 = edge_points_ind + 1 : size(edge_points,2)
        direction2 = edge_points(4,edge_points_ind2); % direction of the second pixel

        % changing direction2 in order to compare
        if direction2 > 0    
            direction2 = direction2 - 180;
        else
            direction2 = direction2 + 180;
        end

        % defining the direction of the line between tow pixels
        m = (edge_points(2,edge_points_ind) - edge_points(2,edge_points_ind2)) / (edge_points(1,edge_points_ind) - edge_points(1,edge_points_ind2));
        m = atand(m);
        m = abs(m);
        n = abs(direction2);
        if n > 90
            n = n - 90;
        end

        % finding pixels on the same circle
        if ((m < (n + degree_tolerance)) && (m > (n - degree_tolerance)))
            % comparing directions
            if ((direction1 < (direction2 + degree_tolerance)) && (direction1 > (direction2 - degree_tolerance)))
                % finding center
                sum_x = edge_points(1,edge_points_ind) + edge_points(1,edge_points_ind2); 
                sum_y = edge_points(2,edge_points_ind) + edge_points(2,edge_points_ind2);
                center(1) = round(sum_x/2);
                center(2) = round(sum_y/2);
                % finding radius
                square_x = (edge_points(1,edge_points_ind) - edge_points(1,edge_points_ind2))^2;
                square_y = (edge_points(2,edge_points_ind) - edge_points(2,edge_points_ind2))^2;
                radius = (square_x + square_y)^0.5 / 2;
                % counting repeats of specific cenetr
                selection_counter(center(1),center(2)) = selection_counter(center(1),center(2)) + 1;
                % inserting amount to the output
                center_radius(center(1),center(2),selection_counter(center(1),center(2))) = radius;
            end
        end
    end                                                                                      
end
close(h)                                                                               %%%%%%%
                                                                                       %%%%%%%
                                                                                       %%%%%%%
%%%%%%%                                       /\                                       %%%%%%%
%%%%%%%                                      /  \                                      %%%%%%%
%%%%%%%                                       ||                                       %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

