function [phi] = a11e02pyramid(i,j,h,order)
    num_elements = int64(2/(h^2));
    num_elements_in_row = 2/h;
    xi = i*h;
    yi = j*h;
    num_points_row = (1/h)+1; 
    element_equations = sym(zeros(num_elements, 1));
    starting_el_index = (j-2)*num_elements+2*(i-2)+1;
    %if selected node is on the boundary
    syms x y
        f1 = 1+(x-i)/h;
        f2 = 1+(y-j)/h;
        f3 = 1-(x-i)/h+(y-j)/h;
        f4 = 1+(x-i)/h-(y-j)/h;
        f5 = 1-(y-j)/h;
        f6 = 1-(x-i)/h;
    if(order > 0)
        f1 = gradient(f1,[x,y]);
        f2 = gradient(f2,[x,y]);
        f3 = gradient(f3,[x,y]);
        f4 = gradient(f4,[x,y]);
        f5 = gradient(f5,[x,y]);
        f6 = gradient(f6,[x,y]);
    end
    %if affected area is on the boundary 
    if ((i == 1 || j == 1) || (i == num_points_row || j == num_points_row) )
       %handle 4 corners 
       if(i == 1 && j == 1)
          element_equations(1,1) = f5;
          element_equations(2,1) = f6;
       elseif(i == 1 && j == num_points_row)
          element_equations((j-2)*num_elements_in_row+1,1) = f3;
       elseif(i == num_points_row && j == num_points_row)
          element_equations((j-1)*num_elements_in_row-1,1) = f1;
          element_equations((j-1)*num_elements_in_row,1) = f2;
       elseif(i == num_points_row && j == 1)
          element_equations(num_elements_in_row,1) = f4;
       end
       %handle 4 situations on the middle of edges
       if((j > 1 && j < num_points_row) && i == 1)
           element_equations((j-2)*num_elements_in_row+1,1) = f3;
           element_equations((j-1)*num_elements_in_row+1,1) = f5;
           element_equations((j-1)*num_elements_in_row+2,1) = f6;     
       elseif((i > 1 && i < num_points_row) && j == 1)
           element_equations((i-1)*2,1) = f4;
           element_equations((i-1)*2+1,1) = f5;
           element_equations((i-1)*2+2,1) = f6;
       elseif((i > 1 && i < num_points_row) && j == num_points_row)
           element_equations((j-2)*num_elements_in_row+2*(i-1)-1,1) = f1;
           element_equations((j-2)*num_elements_in_row+2*(i-1),1) = f2;
           element_equations((j-2)*num_elements_in_row+2*(i-1)+1,1) = f3;
       elseif((j > 1 && j < num_points_row) && i == num_points_row)
           element_equations((j-1)*num_elements_in_row-1,1) = f1;
           element_equations((j-1)*num_elements_in_row,1) = f2;
           element_equations(j*num_elements_in_row,1) = f4;
       end
    else
            element_equations(starting_el_index,1) = f1;
            element_equations(starting_el_index+1,1) = f2;
            element_equations(starting_el_index+2,1) = f3;
            element_equations(starting_el_index+2+(num_elements_in_row-1),1) = f4;
            element_equations(starting_el_index+2+(num_elements_in_row-1)+1,1) = f5;
            element_equations(starting_el_index+2+(num_elements_in_row-1)+2,1) = f6;
    end
     phi = element_equations;
end